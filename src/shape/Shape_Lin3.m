%% Shape_Lin3 Class
%
%% Description
%
% This is a sub-class of the <shape.html Shape> class for the
% implementation of the *Quadratic Line Shape*
% (3-noded isoparametric line).
%
% <<shape_line3.png>>
%
classdef Shape_Lin3 < Shape
    %% Constructor method
    methods
        function this = Shape_Lin3(nodes)
            this = this@Shape(Shape.LIN3,1,3);
            
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
            
            % Vector of local node ids
            this.nodeIdsLcl = [ 1; 3; 2 ];
            
            % Vector of global node ids
            this.nodeIdsGbl = [ nodes(1).id;
                                nodes(3).id;
                                nodes(2).id ];
            
            % Matrix of cartesian coordinates of nodes [X Y Z]
            this.carCoord = [ nodes(1).coord(1) nodes(1).coord(2) nodes(1).coord(3);
                              nodes(2).coord(1) nodes(2).coord(2) nodes(2).coord(3);
                              nodes(3).coord(1) nodes(3).coord(2) nodes(3).coord(3) ];
            
            % Matrix of parametric coordinates of nodes [r]
            this.parCoord = [ -1.0;
                               1.0
                               0.0 ];
            
            % Characteristic length
            p1 = this.carCoord(1,:);
            p2 = this.carCoord(2,:);
            this.len = norm(p1-p2);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtx(this,r,~,~)
            N = this.ShpFcnMtxEdge(r);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxEdge(this,r) 
            N = zeros(1,this.nen);
            
            N(3) = 1-r*r;
            N(1) = 0.5*(1-r) - 0.5*N(3);
            N(2) = 0.5*(1+r) - 0.5*N(3);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxFace(~,~,~)
            N = [];
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxVol(~,~,~,~)
            N = [];
        end
        
        %------------------------------------------------------------------
        function M = MapFcnMtx(this,r,s,t)
            % Isoparametric element shape
            M = this.ShpFcnMtx(r,s,t);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtx(this,r,~,~)
            GradNpar = this.gradShpFcnMtxEdge(r);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxEdge(this,r)
            GradNpar = zeros(this.dim,this.nen);
            
            GradNpar(1,1) = -0.5 + r;
            GradNpar(1,2) =  0.5 + r;
            GradNpar(1,3) = -2*r;
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxFace(~,~,~)
            GradNpar = [];
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxVol(~,~,~,~)
            GradNpar = [];
        end
        
        %------------------------------------------------------------------
        function GradMpar = gradMapFcnMtx(this,r,s,t)
            % Isoparametric element shape
            GradMpar = this.gradShpFcnMtx(r,s,t);
        end
        
        %------------------------------------------------------------------
        function GradMpar = gradMapFcnMtxEdge(this,r)
            % Isoparametric element shape
            GradMpar = this.gradShpFcnMtxEdge(r);
        end
        
        %------------------------------------------------------------------
        function GradMpar = gradMapFcnMtxFace(this,r,s)
            % Isoparametric element shape
            GradMpar = this.gradShpFcnMtxFace(r,s);
        end
        
        %------------------------------------------------------------------
        function GradMpar = gradMapFcnMtxVol(this,r,s,t)
            % Isoparametric element shape
            GradMpar = this.gradShpFcnMtxVol(r,s,t);
        end
        
        
    end
end