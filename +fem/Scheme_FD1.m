%% Scheme_FD1 Class (First-Order Finite Difference Scheme)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Scheme: scheme super-class> to deal
% with first-order finite difference time integration schemes.
%
% The following schemes are implemented with the same expression,
% in which the value of the parameter 'theta' defines each one:
%  * Foward Euler (fully explicit):   theta = 1.0
%  * Backward Euler (fully implicit): theta = 0.0
%  * Crank Nicolson (semi-implicit):  theta = 0.5
%
% Ref.: Fundamentals of the Finite Element Method for Heat and Fluid Flow,
%       R. W. Lewis, P. Nithiarasu and K. N. Seetharamu
%
% Transient semi-discretized system of ODEs on time:
% [C]{Ut} + [K]{U} = {F}
%
% Where:
% {U},{Ut} -> Vector of state variables and their first time derivatives
% [K],[C]  -> Stiffness and "velocity" matrices (constants)
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
%% Class definition
%
classdef Scheme_FD1 < fem.Scheme
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        theta double = double.empty;  % parameter to define scheme
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Scheme_FD1(theta)
            switch theta
                case 0.0
                    type = fem.Scheme.FOWARD_EULER;
                case 1.0
                    type = fem.Scheme.BACKWARD_EULER;
                case 0.5
                    type = fem.Scheme.CRANK_NICOLSON;
            end
            this = this@fem.Scheme(type);
            this.theta = theta;
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Scheme
    methods
        %------------------------------------------------------------------
        % Execute time integration scheme to compute transient solution.
        % Input:
        %  dt:   time step increment (assumed constant)
        %  endt: end time
        %  row:  number of rows for computed arrays (number of free d.o.f.'s)
        %  col:  number of columns for computed arrays (number of time steps)
        %  IC:   matrix of initial conditions (values of state variable only)
        %  K:    global stiffness matrix (free d.o.f.'s only)
        %  C:    global "velocity" matrix (free d.o.f.'s only)
        %  M:    global "acceleration" matrix (free d.o.f.'s only) - not used
        % Output:
        %  U:     matrix of state variable vector for each time step (free d.o.f.'s only)
        %  Ut:    matrix of first time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  Utt:   matrix of second time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  steps: number of computed time steps
        %  times: vector of time values
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
            Y = C - (1-this.theta) * dt * K;
            Z = this.theta*F + (1-this.theta)*F;
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
            
            % Loop through all steps
            for i = dt:dt:endt
                steps = steps + 1;
                times(steps+1) = i;
                
                % Compute state variable at new time step
                U(:,steps+1) = X * (Y * U(:,steps) + dt * Z);
                
                % Compute rate of change at new time step
                Ut(:,steps+1) = Z \ (F - K * U(:,steps+1));
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