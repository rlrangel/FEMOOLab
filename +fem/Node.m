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
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        id    int32                     = int32.empty;                     % identification number
        coord double                    = double.empty;                    % vector of coordinates in global system [X Y Z]
        
        % Boundary conditions
        fixdDOF   logical = logical.empty;      % vector of flags for fixed d.o.f.'s (indicate essential boundary conditions)
        prscDOF   double  = double.empty;       % vector of prescribed values of fixed d.o.f.'s (essential boundary conditions)
        prscNBC   double  = double.empty;       % vector of prescribed values of natural boundary conditions (applied forcing components)
        
        % Initial conditions
        initDOF   double  = double.empty;       % vector of initial values of free d.o.f.'s
        initDOFt  double  = double.empty;       % vector of initial values of 1st time derivatives of free d.o.f.'s
        initDOFtt double  = double.empty;       % vector of initial values of 2nd time derivatives of free d.o.f.'s
        
        % Convection conditions
        convVel   double  = double.empty;       % vector of convection velocity values [vel_X vel_Y vel_Z]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Node()
            return;
        end
    end
end