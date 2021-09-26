%% Material Class
%
%% Description
%
% This is a handle class for the definition of materials.
%
% It stores mechanical and thermal properties of materials.
% All materials are considered to have a linear-elastic bahavior.
% In adition, homogeneous and isotropic properties are considered.
%
classdef Material < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Indentification
        name string = string.empty;   % identification name of material
        
        % Properties
        density   double = double.empty;   % density
        young     double = double.empty;   % Young modulus
        poisson   double = double.empty;   % Poisson ratio
        conduct   double = double.empty;   % thermal conductivity
        hcapacity double = double.empty;   % heat capacity
    end
    
    %% Constructor method
    methods
        function this = Material()
            
        end
    end
end