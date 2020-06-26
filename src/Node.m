%% Node Class
%
% This class defines a node object in the StAnOOP program.
% A node is always considered as a three-dimensional entity, with
% coordinates and boundary conditions specified accordingly.
%
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id         = 0;   % identification number
        coord      = [];  % vector of coordinates in global system [X Y Z]
        ebc        = [];  % vector of essential boundary conditions flags [DX DY DZ]
        prescDispl = [];  % vector of prescribed displacement values [DX DY DZ]
        load       = [];  % vector of applied load components [FX FY FZ]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = Node()
            node.ebc = zeros(3,1);
        end
    end
end