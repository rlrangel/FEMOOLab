%% Material Class
%
%% Description
%
% This class defines a material object of an <element.html element> in the
% FEMOOLab program.
% All materials are considered to have a linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in all
% directions.
%
classdef Material < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id    int32  = int32.empty;     % identification number
        E     double = double.empty;    % elasticity modulus
        v     double = double.empty;    % Poisson ratio
        
        % Thermal properties
        rho   double = double.empty;    % density (specific weight)
        k     double = double.empty;    % thermal conductivity
        cp    double = double.empty;    % specific heat capacity
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function mat = Material()
            return;
        end
    end
end