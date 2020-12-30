%% Node Class
%
%% Description
%
% This class defines a node object in the FEMOOLab program.
% An object of <model.html Model class> has a vector of Node objects.
% Each <shape.html Shape> object of an <edge.html Edge> or <face.html Face>
% object has a Node incidence vector.
% The Node object also points to its connected edges and faces.
%
% A node is always considered as a three-dimensional entity, with
% coordinates and boundary conditions specified accordingly.
% For structural problems, the boundary condition properties are 6-row
% vectors accounting for all displacements/rotations d.o.f.'s:
% [DX DY DZ RX RY RZ]
% For thermal problems, these properties are scalars.
%
%% Class definition
%
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        id    int32    = int32.empty;       % identification number
        ndof  int32    = int32.empty;       % number of total d.o.f.'s (structural = 6; thermal = 1)
        coord double   = double.empty;      % vector of coordinates in global system [X Y Z]
        
        % Incidence
        edges fem.Edge = fem.Edge.empty;    % vector of objects of Edge class (incident edges)
        faces fem.Face = fem.Face.empty;    % vector of objects of Face class (incident faces)
        
        % Boundary conditions
        fixdDOF   logical = logical.empty;  % vector of flags for fixed d.o.f.'s (indicate essential boundary conditions)
        prscDOF   double  = double.empty;   % vector of prescribed values of fixed d.o.f.'s
        initDOF   double  = double.empty;   % vector of initial values of free d.o.f.'s
        initDOFt  double  = double.empty;   % vector of initial values of 1st time derivatives of free d.o.f.'s
        initDOFtt double  = double.empty;   % vector of initial values of 2nd time derivatives of free d.o.f.'s
        prscNBC   double  = double.empty;   % vector of prescribed values of natural boundary conditions (applied forcing components)
        
        % Structural properties (6-row arrays accounting for all d.o.f.'s: [DX DY DZ RX RY RZ])
%         fixDispl logical = logical.empty;  % vector of essential boundary conditions flags
%         ebcDispl double  = double.empty;   % vector of prescribed displacement values
%         iniDispl double  = double.empty;   % vector of initial displacement values
%         iniVeloc double  = double.empty;   % vector of initial velocity values
%         load     double  = double.empty;   % vector of applied load components
        
        % Thermal properties (scalars)
%         fixTemp  logical = logical.empty;  % thermal essential boundary condition flag
%         ebcTemp  double  = double.empty;   % prescribed temperature
%         iniTemp  double  = double.empty;   % initial temperature
%         flux     double  = double.empty;   % applied flux
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Node(ndof)
            this.ndof = ndof;
            this.fixdDOF = false(ndof,1); % By default, all d.o.f.'s are free
        end
    end
end