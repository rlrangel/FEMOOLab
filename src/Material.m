%% Material Class
%
% This class defines a material object in the StAnOOP program.
% All materials are considered to have a linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in all
% directions.
%
classdef Material < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id = 0;  % identification number
        E  = 0;  % elasticity modulus
        v  = 0;  % poisson ratio
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function mat = Material()
            return;
        end
    end
end