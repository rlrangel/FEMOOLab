%% Node Class
%
%% Description
%
% This class defines a node object in the FEMOOLab program.
% An object of <model.html Model class> has a vector of Node objects.
% Each <shape.html Shape> object of an <element.html Element> object
% has a Node incidence vector.
% The Node object also points to its connected elements.
%
% A node is always considered as a three-dimensional entity, with
% coordinates and boundary conditions specified accordingly.
% For structural problems, the boundary condition properties are 6-row
% vectors accounting for all displacements/rotations d.o.f.'s:
% [DX DY DZ RX RY RZ]
% For thermal problems, these properties are scalars (temperature d.o.f.).
%
%% Class definition
%
classdef ExtNode < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        id    int32                     = int32.empty;  % identification number
        coord double                    = double.empty; % vector of coordinates in global system [X Y Z]
        elems                           = [];           % Incident elements: vector of objects of Element_Isoparametric class or
                                                        % Element_Isogeometric class
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = ExtNode()
            return
        end
    end
end