%% Scheme Class
%
%% Description
%
% This is a handle super-class for the definition of numerical schemes for
% time integration.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <Scheme_FD1.html Scheme_FD1>
% * <Scheme_FD2.html Scheme_FD2>
% * <Scheme_RungeKutta.html Scheme_RungeKutta>
% * <Scheme_Newmark.html Scheme_Newmark>
%
classdef Scheme < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of schemes
        FOWARD_EULER       = int32(1);
        BACKWARD_EULER     = int32(2);
        CRANK_NICOLSON     = int32(3);
        FOWARD_GALERKIN    = int32(4);
        BACKWARD_GALERKIN  = int32(5);
        CENTRAL_DIFFERENCE = int32(6);
        RUNGE_KUTTA_1      = int32(7);
        RUNGE_KUTTA_2      = int32(8);
        RUNGE_KUTTA_3      = int32(9);
        RUNGE_KUTTA_4      = int32(10);
        RUNGE_KUTTA_5      = int32(11);
        NEWMARK            = int32(12);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type  int32 = int32.empty;  % flag for type of scheme
        order int32 = int32.empty;  % order of scheme
    end
    
    %% Constructor method
    methods
        function this = Scheme(type,order)
            if (nargin > 0)
                this.type = type;
                this.order = order;
            end            
        end
    end
    
    %% Abstract methods
    methods (Abstract)
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
        [U,Ut,Utt,steps,times] = Execute(this,dt,endt,row,col,IC,K,C,M,F);
    end
end