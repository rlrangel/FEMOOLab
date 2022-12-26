%% Material Class
%
%% Description
%
% This is a handle class for the definition of materials.
%
% It stores thermomechanical properties of materials.
%
classdef Material < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Identification
        name string = string.empty;  % identification name
        
        % Properties
        density      double = double.empty;  % density
        young        double = double.empty;  % Young modulus
        poisson      double = double.empty;  % Poisson ratio
        conductivity double = double.empty;  % thermal conductivity
        hcapacity    double = double.empty;  % heat capacity
    end
    
    %% Constructor method
    methods
        function this = Material()
            return;
        end
    end
end