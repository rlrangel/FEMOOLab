%% Shape Class
%
%% Description
%
% This is a handle super-class for the definition of finite element shapes.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <Shape_Lin2.html Shape_Lin2>
% * <Shape_Lin3.html Shape_Lin3>
% * <Shape_Tri3.html Shape_Tri3>
% * <Shape_Tri6.html Shape_Tri6>
% * <Shape_Quad4.html Shape_Quad4>
% * <Shape_Quad8.html Shape_Quad8>
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
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Identification
        type int32 = int32.empty;  % flag for type of shape
        
        % Dimensions
        dim double = double.empty;  % spatial dimension
        len double = double.empty;  % characteristic length
        
        % Nodes
        n_nodes      int32 = int32.empty;  % number of nodes
        nodes        Node  = Node.empty;   % handles to objects of Node class
        node_ids_lcl int32 = int32.empty;  % vector of local node IDs
        node_ids_gbl int32 = int32.empty;  % vector of global node IDs
        
        % Coordinates
        coord_car double = double.empty;  % matrix of cartesian nodal coordinates
        coord_par double = double.empty;  % matrix of parametric nodal coordinates
        
        % Sub-components
        edge Shape = Shape.empty;  % handle to object of Shape class
    end
    
    %% Constructor method
    methods
        function this = Shape(type,dim,n_nodes)
            if (nargin > 0)
                this.type    = type;
                this.dim     = dim;
                this.n_nodes = n_nodes;
            end
        end
    end
    
    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        initialize(this,nodes);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given element
        % position in parametric coordinates.
        N = ShpFcnMtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given edge
        % position in parametric coordinates.
        N = ShpFcnMtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of geometry map functions at a given element
        % position in parametric coordinates.
        M = MapFcnMtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of d.o.f. shape functions
        % w.r.t. parametric coordinates at a given element position.
        GradNpar = GradShpFcnMtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of d.o.f. shape functions
        % w.r.t. parametric coordinates at a given edge position.
        GradNpar = GradShpFcnMtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of geometry map functions
        % w.r.t. parametric coordinates at a given element position.
        GradMpar = GradMapFcnMtx(this,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of derivatives of geometry map functions
        % w.r.t. parametric coordinates at a given edge position.
        GradMpar = GradMapFcnMtxEdge(this,r);
        
        %------------------------------------------------------------------
        % Returns the local IDs of the nodes of a shape edge for a given
        % pair of corner nodes global IDs.
        % Output:
        %  valid: In case the given corner nodes global IDs do not 
        %         correspond to a valid element shape edge, the returned
        %         value is false. Otherwise, the returned value is true.
        %  n1,n2: The local IDs of the two corner nodes.
        %  mid:   The local ID of the mid node. In case there is no mid
        %         node, it returns zero.
        [valid,n1,n2,mid] = EdgeLocalIds(this,corner1,corner2);
    end
end