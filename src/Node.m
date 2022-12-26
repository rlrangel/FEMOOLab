%% Node Class
%
%% Description
%
% This is a handle class for the definition of nodes.
%
% A node is always considered as a three-dimensional entity,
% with coordinates and conditions specified accordingly.
% For mechanical problems, conditions are 6-row vectors accounting for all
% displacements/rotations d.o.f.'s: [DX DY DZ RX RY RZ].
% For thermal problems, conditions are scalars (temperature d.o.f.).
%
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        id    int32       = int32.empty;        % identification number
        coord double      = double.empty;       % vector of coordinates in global system [X Y Z]
        elems fem.Element = fem.Element.empty;  % vector of objects of Element class (incident elements)
        
        % Initial conditions
        init_dof    double = double.empty;  % vector of initial values of free d.o.f.'s
        init_dof_t  double = double.empty;  % vector of initial values of 1st time derivatives of free d.o.f.'s
        init_dof_tt double = double.empty;  % vector of initial values of 2nd time derivatives of free d.o.f.'s

        % Boundary conditions
        fixed_dof logical = logical.empty;  % vector of flags for fixed d.o.f.'s (indicate essential boundary conditions)
        ebc       double  = double.empty;   % vector of essential boundary conditions (prescribed values of fixed d.o.f.'s)
        nbc       double  = double.empty;   % vector of natural boundary conditions (prescribed values of applied forcing components)
        
        % Convection conditions
        conv_vel double = double.empty;  % vector of convection velocity values [vel_X vel_Y vel_Z]
    end
    
    %% Constructor method
    methods
        function this = Node()
            return;
        end
    end
end