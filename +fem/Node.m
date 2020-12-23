%% Node Class
%
%% Description
%
% This class defines a node object in the StAnOOP program.
% An object of <model.html Model class> has a vector of Node objects.
% Each <shape.html Shape> object of an <element.html Element> object
% has a Node incidence vector.
% A node is always considered as a three-dimensional entity, with
% coordinates and boundary conditions specified accordingly.
%
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id         int32  = int32.empty;   % identification number
        coord      double = double.empty;  % vector of coordinates in global system [X Y Z]
        ebc        int32  = int32.empty;   % vector of essential boundary conditions flags [DX DY DZ RX RY RZ]
        prescDispl double = double.empty;  % vector of prescribed displacement values [DX DY DZ RX RY RZ]
        load       double = double.empty;  % vector of applied load components [FX FY FZ MX MY MZ]
        
        % Thermal properties
        ebc_thermal int32  = int32.empty;   % thermal essential boundary condition flag
        prescTemp   double = double.empty;  % prescribed temperature value
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = Node()
            node.ebc = zeros(6,1);
            node.ebc_thermal = 0;
        end
    end
end