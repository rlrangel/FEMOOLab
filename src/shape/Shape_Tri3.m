%% Shape_Tri3 Class
%
%% Description
%
% This is a sub-class of the <Shape.html Shape> class for the
% implementation of *Linear Triangular Shape*
% (3-noded isoparametric triangle).
%
% <<../images/tutorials/shape_triangle3.png>>
%
classdef Shape_Tri3 < Shape
    %% Constructor method
    methods
        function this = Shape_Tri3(nodes)
            this = this@Shape(Shape.TRI3,2,3);
            
            if (nargin > 0)
                this.initialize(nodes);
            end
        end
    end
    
    %% Public methods: implementation of super-class declarations
    methods
        %------------------------------------------------------------------
        function initialize(this,nodes)
            this.nodes = nodes;
            
            % Vector of local node IDs (ccw order)
            this.node_ids_lcl = [ 1; 2; 3 ];
            
            % Vector of global node IDs (ccw order)
            this.node_ids_gbl = [ nodes(1).id;
                                  nodes(2).id;
                                  nodes(3).id ];
            
            % Matrix of cartesian nodal coordinates [X Y]
            this.coord_car = [ nodes(1).coord(1) nodes(1).coord(2);
                               nodes(2).coord(1) nodes(2).coord(2);
                               nodes(3).coord(1) nodes(3).coord(2) ];
            
            % Matrix of parametric nodal coordinates [r s]
            this.coord_par = [ 0.0  0.0;
                               1.0  0.0;
                               0.0  1.0 ];
            
            % Characteristic length
            x = this.coord_car(:,1);
            y = this.coord_car(:,2);
            area = polyarea(x,y);
            this.len = area^(1/2);

            % Sub-components
            this.edge = Shape_Lin2();
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtx(this,r,s)
            N = zeros(1,this.n_nodes);
            
            N(1) = 1 - r - s;
            N(2) = r;
            N(3) = s;
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxEdge(this,r) 
            N = this.edge.ShpFcnMtx(r,[]);
        end
        
        %------------------------------------------------------------------
        function M = MapFcnMtx(this,r,s)
            % Isoparametric element shape
            M = this.ShpFcnMtx(r,s);
        end
        
        %------------------------------------------------------------------
        function GradNpar = GradShpFcnMtx(this,~,~)
            GradNpar = zeros(this.dim,this.n_nodes);
            
            GradNpar(1,1) = -1.0;
            GradNpar(2,1) = -1.0;
            GradNpar(1,2) =  1.0;
            GradNpar(2,2) =  0.0;
            GradNpar(1,3) =  0.0;
            GradNpar(2,3) =  1.0;
        end
        
        %------------------------------------------------------------------
        function GradNpar = GradShpFcnMtxEdge(this,r)
             GradNpar = this.edge.GradShpFcnMtx(r,[]);
        end
        
        %------------------------------------------------------------------
        function GradMpar = GradMapFcnMtx(this,r,s)
            % Isoparametric element shape
            GradMpar = this.GradShpFcnMtx(r,s);
        end
        
        %------------------------------------------------------------------
        function GradMpar = GradMapFcnMtxEdge(this,r)
            % Isoparametric element shape
            GradMpar = this.GradShpFcnMtxEdge(r);
        end
        
        %------------------------------------------------------------------
        function [valid,n1,n2,mid] = EdgeLocalIds(this,corner1,corner2)
            valid = false;
            n1    = 0;
            n2    = 0;
            mid   = 0;

            % Get IDs of corner nodes and check for consistency
            for i = 1:this.n_nodes
                node = this.nodes(i).id;
                if node == corner1
                    n1 = i;
                    break;
                end
            end
            
            if n1 == 0
                return;
            end
            
            for i = 1:this.n_nodes
                node = this.nodes(i).id;
                if node == corner2
                    n2 = i;
                    break;
                end
            end
            
            if n2 == 0
                n1 = 0;
                return;
            end
            
            valid = true;
        end
    end
end