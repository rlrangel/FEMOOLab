%% Scheme_FD2 Class (Second-Order Finite Difference Scheme)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Scheme: scheme super-class> to deal
% with second-order finite difference time integration schemes.
%
% The second-order central difference explicit scheme is implemented for
% solving a system of ODEs of second-order.
%
% Ref.: A survey of direct time-integration methods in computational
%       structural dynamics - I: Explicit methods,
%       M. A. Dokainish and K. Subbara
%       Computers & Structures Vol. 32., No. 6., pp. 1371-1386, 1989
%
% Transient semi-discretized system of ODEs on time:
% [M]{Utt} + [C]{Ut} + [K]{U} = {F}
%
% Where:
% {U},{Ut},{Utt} -> Vector of state variables and their 1st and 2nd time derivatives
% [K],[C],[M]    -> corresponding matrices (constants)
% {F}            -> Forcing vector (currently assumed constant on time)
%
% Fully discretized system using second-order central finite difference
% approximation on time:
% [X]{U}i+1 = {F} - [Y]{U}i - [Z]{U}i-1
%
% Where:
% [X] = 1/(dt^2)*[M] + 1/(2dt)*[C]
% [Y] = K - 2/(dt^2)*[M]
% [Z] = 1/(dt^2)*[M] - 1/(2dt)*[C]
%
%% Class definition
%
classdef Scheme_FD2 < fem.Scheme
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Scheme_FD2()
            this = this@fem.Scheme(fem.Scheme.CENTRAL_DIFFERENCE,2);
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
        %  col:  number of columns for computed arrays (number of time steps + init. cond.)
        %  IC:   matrix of initial conditions
        %  K:    global stiffness matrix (free d.o.f.'s only)
        %  C:    global matrix related to 1st time derivative (free d.o.f.'s only)
        %  M:    global matrix related to 2nd time derivative (free d.o.f.'s only)
        %  F:    global forcing vector (free d.o.f.'s only) - currently assumed constant!
        % Output:
        %  U:     matrix of state variable vector for each time step (free d.o.f.'s only)
        %  Ut:    matrix of first time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  Utt:   matrix of second time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  steps: number of computed time steps
        %  times: vector of time values
        function [U,Ut,Utt,steps,times] = execute(~,dt,endt,row,col,IC,K,C,M,F)
            % Initialize solution matrices
            U   = zeros(row,col);
            Ut  = zeros(row,col);
            Utt = zeros(row,col);
            
            % Set initial condition of state variable and first time derivative
            U(:,1)  = IC(:,1);
            Ut(:,1) = IC(:,2);
            
            % Compute initial condition of second time derivative from global system:
            % [M]*{Utt} = {F} - [K]{U} - [C]*{Ut}
            Utt(:,1) = M \ (F - K * U(:,1) - C * Ut(:,1));
            
            % Compute auxiliar arrays
            X = inv(1/(dt^2) * M + 1/(2*dt) * C);
            Y = K - 2/(dt^2) * M;
            Z = 1/(dt^2) * M - 1/(2*dt) * C;
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
            
            % Loop through all steps
            for i = dt:dt:double(endt)
                steps = steps + 1;
                times(steps+1) = i;
                
                % Get state variable at previous time step
                u = U(:,steps) - dt * Ut(:,steps) + dt^2/2 * Utt(:,steps);
                
                % Compute state variable at new time step
                U(:,steps+1) = X * (F - Y * U(:,steps) - Z * u);
                
                % Compute first and second derivatives at new time step
                Ut(:,steps+1)  = (U(:,steps+1) - u)/(2*dt);
                Utt(:,steps+1) = (U(:,steps+1) - 2*U(:,steps) + u)/(dt^2);
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