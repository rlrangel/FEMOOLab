%% Node class
classdef Node < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 0;                    % identification number
        coords = [];               % 3 position vector of nodal coordinates [X, Y, Z]
        ebc = zeros(1,6);          % 6 position vector of essential boundary condition flags [0=free, 1=fixed]
        nodalLoad = zeros(1,6);    % 6 position vector of nodal load components [fx, fy, fz, mx, my, mz]
        prescDispl = zeros(1,6);   % 6 position vector of prescribed displacements [dx, dy, dz, rx, ry, rz]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = Node(id,coords,ebc,nl,presdispl)
            if (nargin > 0)
                node.id = id;
                node.coords = coords;
                node.ebc = ebc;
                node.nodalLoad = nl;
                node.prescDispl = presdispl;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a Node object.
        function clean(node)
            node.id = 0;
            node.coords = [];
            node.ebc = [];
            node.nodalLoad = [];
            node.prescDispl = [];
        end
    end
end