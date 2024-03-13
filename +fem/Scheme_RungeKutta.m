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
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Scheme_RungeKutta(order)
            switch order
                case 1
                    type = fem.Scheme.RUNGE_KUTTA_1;
                case 2
                    type = fem.Scheme.RUNGE_KUTTA_2;
                case 3
                    type = fem.Scheme.RUNGE_KUTTA_3;
                case 4
                    type = fem.Scheme.RUNGE_KUTTA_4;
                case 5
                    type = fem.Scheme.RUNGE_KUTTA_5;
            end
            this = this@fem.Scheme(type,order);
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
        function [U,Ut,Utt,steps,times] = execute(this,dt,endt,row,col,IC,K,C,M,F)
            % Call method for corresponding order
            switch this.order
                case 1
                    [U,Ut,Utt,steps,times] = this.RK1(dt,endt,row,col,IC,K,C,M,F);
                case 2
                    [U,Ut,Utt,steps,times] = this.RK2(dt,endt,row,col,IC,K,C,M,F);
                case 3
                    [U,Ut,Utt,steps,times] = this.RK3(dt,endt,row,col,IC,K,C,M,F);
                case 4
                    [U,Ut,Utt,steps,times] = this.RK4(dt,endt,row,col,IC,K,C,M,F);
                case 5
                    [U,Ut,Utt,steps,times] = this.RK5(dt,endt,row,col,IC,K,C,M,F);
            end
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Fourth-order Runge-Kutta scheme.
        function [U,Ut,Utt,steps,times] = RK1(dt,endt,row,col,IC,K,C,M,F)
            
            % NOT IMPLEMENTED !!!
            
            % Initialize solution matrices
            U   = zeros(row,col);
            Ut  = zeros(row,col);
            Utt = zeros(row,col);
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
        end
        
        %------------------------------------------------------------------
        % Fourth-order Runge-Kutta scheme.
        function [U,Ut,Utt,steps,times] = RK2(dt,endt,row,col,IC,K,C,M,F)
            
            % NOT IMPLEMENTED !!!
            
            % Initialize solution matrices
            U   = zeros(row,col);
            Ut  = zeros(row,col);
            Utt = zeros(row,col);
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
        end
        
        %------------------------------------------------------------------
        % Fourth-order Runge-Kutta scheme.
        function [U,Ut,Utt,steps,times] = RK3(dt,endt,row,col,IC,K,C,M,F)
            
            % NOT IMPLEMENTED !!!
            
            % Initialize solution matrices
            U   = zeros(row,col);
            Ut  = zeros(row,col);
            Utt = zeros(row,col);
            
            % Initialize step and time counting
            steps = 0;
            times = zeros(1,col);
        end
        
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
        
        %------------------------------------------------------------------
        % Fourth-order Runge-Kutta scheme.
        function [U,Ut,Utt,steps,times] = RK5(dt,endt,row,col,IC,K,C,M,F)
            
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