%% Shape Class
%
%% Description
%
% This is an abstract super-class that generically specifies an
% <element.html element> shape in the FEMOOLab program.
% Essentially, this super-class declares abstract methods that define the
% general behavior of elements. These abstract methods are the functions
% that should be implemented in a derived sub-class that deals with
% specific types of elements.
%
%% Subclasses
%
% * <shape_tria3.html Shape_Tria3: 3-node linear triangle shape subclass>
% * <shape_quad4.html Shape_Quad4: 4-node bilinear quadrilateral shape subclass>
% * <shape_tria6.html Shape_Tria6: 6-node quadratic triangle shape subclass>
% * <shape_quad8.html Shape_Quad8: 8-node quadratic quadrilateral shape subclass>
%
%% Class definition
%
classdef Shape < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of shape
        TRIA3 = int32(1);
        QUAD4 = int32(2);
        TRIA6 = int32(3);
        QUAD8 = int32(4);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        type  int32  = int32.empty;                 % flag for type of shape
        order int32  = int32.empty;                 % linear (1) or quadratic (2)
        dim   double = double.empty;                % shape dimension
        size  double = double.empty;                % shape size (length, area, or volume)
        
        % Coordinates
        carCoord        double   = double.empty;    % matrix of cartesian nodal coordinates
        parCoord        double   = double.empty;    % matrix of parametric nodal coordinates
        
        % Nodes
        nen             int32    = int32.empty;     % number of nodes
        nodes           fem.Node = fem.Node.empty;  % vector of objects of Node class
        ccwLocalNodeIds int32    = int32.empty;     % vector of local node ids in ccw order
        ccwNodeIds      int32    = int32.empty;     % vector of global node ids in ccw order
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape(type,dim,nen)
            this.type = type;
            this.dim  = dim;
            this.nen  = nen;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions at a given position in
        % parametric coordinates.
        M = Mmtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given position in
        % parametric coordinates.
        N = Nmtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge d.o.f. shape functions at a given
        % position in parametric coordinates.
        N = NmtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        GradMpar = gradMmtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        GradNpar = gradNmtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        GradMpar = gradMmtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Returns the local ids of the nodes of a shape edge for a given
        % pair of corner nodes global ids.
        % In case the given corner nodes global ids do not correspond to
        % a valid element shape edge, the returned value for the valid
        % parameter is false. Otherwise, the returned value is true.
        % The two corner nodes local ids are returned in parameters n1 and
        % n2. In case there is a mid node at the target edge, the local id
        % of this node is returned in the parameter mid. Otherwise, the
        % returned value of mid is zero.
        [valid,n1,n2,mid] = edgeLocalIds(this,corner1,corner2);
    end
end