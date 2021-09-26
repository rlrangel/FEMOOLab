%% Shape Class
%
%% Description
%
% This is a handle super-class for the definition of finite element shapes.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <shape_lin2.html Shape_Lin2>
% * <shape_lin3.html Shape_Lin3>
% * <shape_tri3.html Shape_Tri3>
% * <shape_tri6.html Shape_Tri6>
% * <shape_quad4.html Shape_Quad4>
% * <shape_quad8.html Shape_Quad8>
% * <shape_tet4.html Shape_Tet4>
% * <shape_tet10.html Shape_Tet10>
% * <shape_hex8.html Shape_Hex8>
% * <shape_hex20.html Shape_Hex20>
%
classdef Shape < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of shape
        LIN2  = int32(1);
        LIN3  = int32(2);
        TRI3  = int32(3);
        TRI6  = int32(4);
        QUAD4 = int32(5);
        QUAD8 = int32(6);
        TET4  = int32(7);
        TET10 = int32(8);
        HEX8  = int32(9);
        HEX20 = int32(10);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Indentification
        type int32 = int32.empty;   % flag for type of shape
        
        % Dimensions
        dim double = double.empty;   % spatial dimension
        len double = double.empty;   % characteristic length
        
        % Nodes
        nen        int32 = int32.empty;   % number of nodes
        nodes      Node  = Node.empty;    % handles to objects of Node class
        nodeIdsLcl int32 = int32.empty;   % vector of local node ids
        nodeIdsGbl int32 = int32.empty;   % vector of global node ids
        
        % Coordinates
        carCoord double = double.empty;   % matrix of cartesian nodal coordinates
        parCoord double = double.empty;   % matrix of parametric nodal coordinates
        
        % Sub-components
        edge Shape = Shape.empty;   % handle to object of Shape class
        face Shape = Shape.empty;   % handle to object of Shape class
    end
    
    %% Constructor method
    methods
        function this = Shape(type,dim,nen)
            if (nargin > 0)
                this.type = type;
                this.dim  = dim;
                this.nen  = nen;
            end
        end
    end
    
    %% Abstract methods: implemented in derived sub-classes
    methods (Abstract)
        %------------------------------------------------------------------
        initialize(this,nodes);
        
        %------------------------------------------------------------------
        % Evaluate matrix of element d.o.f. shape functions at a given
        % position in parametric coordinates.
        N = ShpFcnMtx(this,r,s,t);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given edge
        % position in parametric coordinates.
        N = ShpFcnMtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given face
        % position in parametric coordinates.
        N = ShpFcnMtxFace(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given volume
        % position in parametric coordinates.
        N = ShpFcnMtxVol(this,r,s,t);
        
        %------------------------------------------------------------------
        % Evaluate matrix of element geometry map functions at a given
        % position in parametric coordinates.
        M = MapFcnMtx(this,r,s,t);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of element d.o.f. shape functions
        % w.r.t. parametric coordinates at a given position.
        GradNpar = gradShpFcnMtx(this,r,s,t);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of d.o.f. shape functions
        % w.r.t. parametric coordinates at a given edge position.
        GradNpar = gradShpFcnMtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of d.o.f. shape functions
        % w.r.t. parametric coordinates at a given face position.
        GradNpar = gradShpFcnMtxFace(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of d.o.f. shape functions
        % w.r.t. parametric coordinates at a given volume position.
        GradNpar = gradShpFcnMtxVol(this,r,s,t);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of element geometry map functions
        % w.r.t. parametric coordinates at a given position.
        GradMpar = gradMapFcnMtx(this,r,s,t);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of geometry map functions
        % w.r.t. parametric coordinates at a given edge position.
        GradMpar = gradMapFcnMtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of geometry map functions
        % w.r.t. parametric coordinates at a given face position.
        GradMpar = gradMapFcnMtxFace(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of geometry map functions
        % w.r.t. parametric coordinates at a given volume position.
        GradMpar = gradMapFcnMtxVol(this,r,s,t);
        
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