%% Scheme_FD1 Class
%
%% Description
%
% This is a sub-class of the <Scheme.html Scheme> class for the
% implementation of *First-Order Finite Difference Scheme* of numerical
% time integration.
%
% The following schemes are implemented with the same expression,
% in which the value of the parameter 'theta' defines each one:
%  * Foward Euler (fully explicit):      theta = 1.0
%  * Backward Euler (fully implicit):    theta = 0.0
%  * Crank Nicolson (semi-implicit):     theta = 1/2
%  * Foward Galerkin:                    theta = 1/3
%  * Backward Galerkin:                  theta = 2/3
%
% Ref.:
% R. W. Lewis, P. Nithiarasu and K. N. Seetharamu -
% Fundamentals of the Finite Element Method for Heat and Fluid Flow.
%
% Transient semi-discretized system of ODEs on time:
% [C]{Ut} + [K]{U} = {F}
%
% Where:
% {U},{Ut} -> Vector of state variables and their first time derivatives
% [K],[C]  -> corresponding matrices (constants)
% {F}      -> Forcing vector (currently assumed constant on time)
%
% Fully discretized system using first-order finite difference
% approximation on time:
% [X]{U}i+1 = [Y]{U}i + dt*{Z}
%
% Where:
% [X] = [C] + theta * dt * [K]
% [Y] = [C] - (1-theta) * dt * [K]
% {Z} = theta*{F}i+1 + (1-theta){F}i
%
classdef Scheme_FD1 < Scheme
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        theta double = double.empty;  % parameter to define scheme
    end
    
    %% Constructor method
    methods
        function this = Scheme_FD1(theta)
            type = Scheme.FOWARD_EULER; % default scheme
            switch theta
                case 0.0
                    type = Scheme.FOWARD_EULER;
                case 1.0
                    type = Scheme.BACKWARD_EULER;
                case 0.5
                    type = Scheme.CRANK_NICOLSON;
                case 1/3
                    type = Scheme.FOWARD_GALERKIN;
                case 2/3
                    type = Scheme.BACKWARD_GALERKIN;
            end
            this = this@Scheme(type,1);
            this.theta = theta;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function [U,Ut,Utt,steps,times] = execute(this,dt,endt,row,col,IC,K,C,~,F)
            % Initialize solution matrices
            U   = zeros(row,col);
            Ut  = zeros(row,col);
            Utt = zeros(row,col); % second time derivative not computed on first order schemes
            
            % Set initial condition of state variable
            U(:,1) = IC(:,1);
            
            % Compute initial condition of rate of change from global system:
            % [C]*{Ut} = {F} - [K]{U}
            Ut(:,1) = C \ (F - K * U(:,1));
            
            % Compute auxiliar arrays
            X = inv(C + this.theta * dt * K);
            Y = C - (1 - this.theta) * dt * K;
            Z = this.theta * F + (1 - this.theta) * F;
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
            
            % Loop through all steps
            for i = dt:dt:double(endt)
                steps = steps + 1;
                times(steps+1) = i;
                
                % Compute state variable at new time step
                U(:,steps+1) = X * (Y * U(:,steps) + dt * Z);
                
                % Compute rate of change at new time step
                Ut(:,steps+1) = C \ (F - K * U(:,steps+1));
            end
            
            % Clean unused steps
            if (steps < col-1)
                U     = U(:,1:steps+1);
                Ut    = Ut(:,1:steps+1);
                Utt   = Utt(:,1:steps+1);
                times = times(1:steps+1);
            end
        end
    end
end