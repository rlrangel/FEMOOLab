%% Scheme_RungeKutta Class
%
%% Description
%
% This is a sub-class of the <Scheme.html Scheme> class for the
% implementation of *Runge Kutta Scheme* of numerical time integration.
%
classdef Scheme_RungeKutta < Scheme
    %% Constructor method
    methods
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
            this = this@Scheme(type,order);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
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
        % First-order Runge-Kutta scheme.
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
        % Second-order Runge-Kutta scheme.
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
        % Third-order Runge-Kutta scheme.
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
        % Fifth-order Runge-Kutta scheme.
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