%% Scheme_Newmark Class
%
%% Description
%
% This is a sub-class of the <Scheme.html Scheme> class for the
% implementation of *Newmark Scheme* of numerical time integration.
%
classdef Scheme_Newmark < Scheme
    %% Constructor method
    methods
        function this = Scheme_Newmark()
            this = this@Scheme(Scheme.NEWMARK,2);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
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