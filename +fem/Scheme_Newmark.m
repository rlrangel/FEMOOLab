%% Scheme_Newmark Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Scheme: scheme super-class> to deal
% with Newmark time integration scheme.
%
% ADD METHODOLOGY AND REFERENCES HERE !!!
%
%% Class definition
%
classdef Scheme_Newmark < fem.Scheme
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Scheme_Newmark()
            this = this@fem.Scheme(fem.Scheme.NEWMARK,2);
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