%% Node_Isogeometric Class
%
%% Description
%
% This class defines a control point object in the FEMOOLab program.
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

classdef Node_Isogeometric_Bezier_Extraction < fem.Node
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        weight double                  = double.empty;                     % vector of weights
        elems fem.Element_Isogeometric_Bezier_Extraction = fem.Element_Isogeometric_Bezier_Extraction.empty; % vector of objects of Element_Isogeometric_Bezier_Extraction class (incident elements)
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Node_Isogeometric()
            return;
        end
    end
end