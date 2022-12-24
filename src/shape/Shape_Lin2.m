%% Shape_Lin2 Class
%
%% Description
%
% This is a sub-class of the <Shape.html Shape> class for the
% implementation of *Linear Line Shape*
% (2-noded isoparametric line).
%
classdef Shape_Lin2 < Shape
    %% Constructor method
    methods
        function this = Shape_Lin2(nodes)
            this = this@Shape(Shape.LIN2,1,2);
            
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
            
            % Vector of local node IDs
            this.node_ids_lcl = [ 1; 2 ];
            
            % Vector of global node IDs
            this.node_ids_gbl = [ nodes(1).id;
                                  nodes(2).id ];
            
            % Matrix of cartesian nodal coordinates [X Y]
            this.coord_car = [ nodes(1).coord(1) nodes(1).coord(2);
                               nodes(2).coord(1) nodes(2).coord(2) ];
            
            % Matrix of parametric nodal coordinates [r]
            this.coord_par = [ -1.0;
                                1.0 ];
            
            % Characteristic length
            p1 = this.coord_car(1,:);
            p2 = this.coord_car(2,:);
            this.len = norm(p1-p2);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtx(this,r,~)
            N = zeros(1,this.n_nodes);
            
            N(1) = 0.5 * (1 - r);
            N(2) = 0.5 * (1 + r);
        end
        
        %------------------------------------------------------------------
        function M = MapFcnMtx(this,r,s)
            % Isoparametric element shape
            M = this.ShpFcnMtx(r,s);
        end
        
        %------------------------------------------------------------------
        function GradNpar = GradShpFcnMtx(this,~,~)
            GradNpar = zeros(this.dim,this.n_nodes);
            
            GradNpar(1,1) = -0.5;
            GradNpar(1,2) =  0.5;
        end
        
        %------------------------------------------------------------------
        function GradMpar = GradMapFcnMtx(this,r,s)
            % Isoparametric element shape
            GradMpar = this.GradShpFcnMtx(r,s);
        end
    end
end