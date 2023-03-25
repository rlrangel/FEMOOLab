%% Anme Class (Analysis Method)
%
%% Description
%
% This is an abstract class that generically specifies an analysis 
% Method in the FEMOOLab program.
% Essentially, the analysis method refers whether the analysis
% is in the isoparametric or isogeometric scope.
%
%% Class definition
%
classdef Anme < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of analysis Method
        ISOPARAMETRIC = int32(1);
        ISOGEOMETRIC  = int32(2);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;  % flag for analysis method
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anme(type)
            this.type = type;
        end
    end
end