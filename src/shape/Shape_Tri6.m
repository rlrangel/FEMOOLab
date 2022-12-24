%% Shape_Tri6 Class
%
%% Description
%
% This is a sub-class of the <Shape.html Shape> class for the
% implementation of *Quadratic Triangular Shape*
% (6-noded isoparametric triangle).
%
% <<../images/tutorials/shape_triangle6.png>>
%
classdef Shape_Tri6 < Shape
    %% Constructor method
    methods
        function this = Shape_Tri6(nodes)
            this = this@Shape(Shape.TRI6,2,6);
            
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
            this.node_ids_lcl = [ 1; 4; 2; 5; 3; 6 ];
            
            % Vector of global node IDs (ccw order)
            this.node_ids_gbl = [ nodes(1).id;
                                  nodes(4).id;
                                  nodes(2).id;
                                  nodes(5).id;
                                  nodes(3).id;
                                  nodes(6).id ];
            
            % Matrix of cartesian nodal coordinates [X Y]
            this.coord_car = [ nodes(1).coord(1) nodes(1).coord(2);
                               nodes(2).coord(1) nodes(2).coord(2);
                               nodes(3).coord(1) nodes(3).coord(2);
                               nodes(4).coord(1) nodes(4).coord(2);
                               nodes(5).coord(1) nodes(5).coord(2);
                               nodes(6).coord(1) nodes(6).coord(2) ];
            
            % Matrix of parametric nodal coordinates [r s]
            this.coord_par = [ 0.0  0.0;
                               1.0  0.0;
                               0.0  1.0;
                               0.5  0.0;
                               0.5  0.5;
                               0.0  0.5 ];
            
            % Characteristic length
            x = this.coord_car(:,1);
            y = this.coord_car(:,2);
            area = polyarea(x,y);
            this.len = area^(1/2);

            % Sub-components
            this.edge = Shape_Lin3();
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtx(this,r,s)
            N = zeros(1,this.n_nodes);
            
            N(1) = 1 - 3*r - 3*s + 4*r*s + 2*r*r + 2*s*s;
            N(2) = 2*r*r - r;
            N(3) = 2*s*s - s;
            N(4) = 4*r - 4*r*r - 4*r*s;
            N(5) = 4*r*s;
            N(6) = 4*s - 4*r*s - 4*s*s;
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
        function GradNpar = GradShpFcnMtx(this,r,s)
            GradNpar = zeros(this.dim,this.n_nodes);
            
            GradNpar(1,1) =  4.0 * r + 4.0 * s - 3.0;
            GradNpar(2,1) =  4.0 * r + 4.0 * s - 3.0;
            GradNpar(1,2) =  4.0 * r - 1.0;
            GradNpar(2,2) =  0.0;
            GradNpar(1,3) =  0.0;
            GradNpar(2,3) =  4.0 * s - 1.0;
            GradNpar(1,4) =  4.0 - 8.0 * r - 4.0 * s;
            GradNpar(2,4) = -4.0 * r;
            GradNpar(1,5) =  4.0 * s;
            GradNpar(2,5) =  4.0 * r;
            GradNpar(1,6) = -4.0 * s;
            GradNpar(2,6) =  4.0 - 4.0 * r - 8.0 * s;
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

            % Get ids of corner nodes and check for consistency
            for i = 1:3
                node = this.nodes(i).id;
                if node == corner1
                    n1 = i;
                    break;
                end
            end
            
            if n1 == 0
                return;
            end
            
            for i = 1:3
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
           
            % Check for corner node consistency and get mid point id
            if n1 == 1
                if n2 == 2
                    mid = 4;
                else % if n2 == 3
                    mid = 6;
                end
            elseif n1 == 2
                if n2 == 1
                    mid = 4;
                else % if n2 == 3
                    mid = 5;
                end
            else % if n1 == 3
                if n2 == 1
                    mid = 6;
                else % if n2 == 2
                    mid = 5;
                end
            end
            
            valid = true;
        end
    end
end