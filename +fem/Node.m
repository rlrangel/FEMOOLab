%% Node Class
%
%% Description
%
% This class defines a node object in the FEMOOLab program.
% An object of <model.html Model class> has a vector of Node objects.
% Each <shape.html Shape> object of an <element.html Element> object
% has a Node incidence vector.
% A node is always considered as a three-dimensional entity, with
% coordinates and boundary conditions specified accordingly.
%
%% Class definition
%
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        id       int32  = int32.empty;     % identification number
        coord    double = double.empty;    % vector of coordinates in global system [X Y Z]
        
        % Structural properties (6-row arrays accounting for all d.o.f.'s: [DX DY DZ RX RY RZ])
        fixDispl logical = logical.empty;  % vector of essential boundary conditions flags 
        ebcDispl double  = double.empty;   % vector of prescribed displacement values
        iniDispl double  = double.empty;   % vector of initial displacement values
        iniVeloc double  = double.empty;   % vector of initial velocity values
        load     double  = double.empty;   % vector of applied load components
        
        % Thermal properties
        fixTemp  logical = logical.empty;  % thermal essential boundary condition flag
        ebcTemp  double  = double.empty;   % prescribed temperature
        iniTemp  double  = double.empty;   % initial temperature
        flux     double  = double.empty;   % applied flux
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Node()
            this.fixDispl = false(6,1);
            this.fixTemp  = false;
        end
    end
end