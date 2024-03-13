%% Material class
classdef Material < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 0;           % identification number
        elasticity = 0;   % elasticity modulus
        poisson = 0;      % poisson ratio
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function material = Material(id,e,v)
            if (nargin > 0)
                material.id = id;
                material.elasticity = e;
                material.poisson = v;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a Material object.
        function clean(material)
            material.id = 0;
            material.elasticity = 0;
            material.poisson = 0;
        end
    end
end