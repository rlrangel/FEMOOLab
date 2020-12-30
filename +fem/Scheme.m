%% Scheme Class
%
%% Description
%
% This is an abstract super-class that generically specifies a time 
% integration scheme in the FEMOOLab program.
% Essentially, this super-class declares abstract methods that define
% the general behavior of time integration scheme types. These abstract
% methods are the functions that should be implemented in a derived
% sub-class that deals with specific types of scheme.
%
%% Subclasses
%
% * <scheme_fd1.html Scheme_FD1: first-order finite difference scheme subclass>
% * <scheme_fd2.html Scheme_FD2: second-order finite difference scheme subclass>
% * <scheme_rungekutta.html Scheme_RungeKutta: Runge-Kutta scheme subclass>
% * <scheme_newmark.html Scheme_Newmark: Newmark scheme subclass>
%
%% Class definition
%
classdef Scheme < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of schemes
        FOWARD_EULER       = int32(1);
        BACKWARD_EULER     = int32(2);
        CRANK_NICOLSON     = int32(3);
        CENTRAL_DIFFERENCE = int32(4);
        RUNGE_KUTTA_1      = int32(5);
        RUNGE_KUTTA_2      = int32(6);
        RUNGE_KUTTA_3      = int32(7);
        RUNGE_KUTTA_4      = int32(8);
        RUNGE_KUTTA_5      = int32(9);
        NEWMARK            = int32(10);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;  % flag for type of scheme
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Scheme(type)
            this.type = type;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
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
        %  C:    global "velocity" matrix (free d.o.f.'s only)
        %  M:    global "acceleration" matrix (free d.o.f.'s only)
        %  F:    global forcing vector (free d.o.f.'s only) - currently assumed constant!
        % Output:
        %  U:     matrix of state variable vector for each time step (free d.o.f.'s only)
        %  Ut:    matrix of first time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  Utt:   matrix of second time derivative of state variable vector for each time step (free d.o.f.'s only)
        %  steps: number of computed time steps
        %  times: vector of time values
        [U,Ut,Utt,steps,times] = execute(this,dt,endt,row,col,IC,K,C,M,F);
    end
end