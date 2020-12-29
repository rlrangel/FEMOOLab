%% Scheme_RungeKutta Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Scheme: scheme super-class> to deal
% with Runge-Kutta time integration schemes.
%
% ADD METHODOLOGY AND REFERENCES HERE !!!
%
%% Class definition
%
classdef Scheme_RungeKutta < fem.Scheme
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        order int32 = int32.empty;  % order of Runge-Kutta scheme
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Scheme_RungeKutta(order)
            switch order
                case 4
                    type = fem.Scheme.RUNGE_KUTTA_4;
            end
            this = this@fem.Scheme(type);
            this.order = order;
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
        %  C:    global "velocity" matrix (free d.o.f.'s only)
        %  M:    global "acceleration" matrix (free d.o.f.'s only)
        % Output:
        %  U:     matrix of state variable vector for each time step (free d.o.f.'s only)
        %  Ut:    matrix of first time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  Utt:   matrix of second time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  steps: number of computed time steps
        %  times: vector of time values
        function [U,Ut,Utt,steps,times] = execute(this,dt,endt,row,col,IC,K,C,M,F)
            % Call method for corresponding order
            if (this.order == 4)
                [U,Ut,Utt,steps,times] = this.RK4(dt,endt,row,col,IC,K,C,M,F);
            end
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Fourth-order Runge-Kutta scheme.
        function [U,Ut,Utt,steps,times] = RK4(dt,endt,row,col,IC,K,C,M,F)
            
            % NOT IMPLEMENTED !!!
            
            % Initialize solution matrices
            U   = zeros(row,col);
            Ut  = zeros(row,col);
            Utt = zeros(row,col);
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
        end
    end
end