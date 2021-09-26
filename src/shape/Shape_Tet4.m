%% Shape_Tet4 Class
%
%% Description
%
% This is a sub-class of the <shape.html Shape> class for the
% implementation of the *Linear Thetahedral Shape*
% (4-noded isoparametric tetrahedron).
%
% <<shape_tetrahedron4.png>>
%
classdef Shape_Tet4 < Shape
    %% Constructor method
    methods
        function this = Shape_Tet4(nodes)
            this = this@Shape(Shape.TET4,3,4);
            
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
            this.nodeIdsLcl = [ 1; 2; 3; 4 ];
            
            % Vector of global node ids
            this.nodeIdsGbl = [ nodes(1).id;
                                nodes(2).id;
                                nodes(3).id;
                                nodes(4).id ];
            
            % Matrix of cartesian coordinates of nodes [X Y Z]
            this.carCoord = [ nodes(1).coord(1) nodes(1).coord(2) nodes(1).coord(3);
                              nodes(2).coord(1) nodes(2).coord(2) nodes(2).coord(3);
                              nodes(3).coord(1) nodes(3).coord(2) nodes(3).coord(3);
                              nodes(4).coord(1) nodes(4).coord(2) nodes(4).coord(3) ];
            
            % Matrix of parametric coordinates of nodes [r s t]
            this.parCoord = [ 0.0  0.0  0.0;
                              1.0  0.0  0.0;
                              0.0  1.0  0.0;
                              0.0  0.0  1.0 ];
            
            % Sub-components
            this.edge = Shape_Lin2();
            this.face = Shape_Tri3();
            
            % Characteristic length
            x = this.carCoord(:,1);
            y = this.carCoord(:,2);
            z = this.carCoord(:,3);
            [~,vol] = convhull(x,y,z);
            this.len = vol^(1/3);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtx(this,r,s,t)
            N = this.ShpFcnMtxVol(r,s,t);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxEdge(this,r) 
            N = this.edge.ShpFcnMtx(r,[],[]);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxFace(this,r,s)
            N = this.face.ShpFcnMtx(r,s,[]);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxVol(this,r,s,t)
            N = zeros(1,this.nen);
            
            N(1) = 1 - r - s - t;
            N(2) = r;
            N(3) = s;
            N(4) = t;
        end
        
        %------------------------------------------------------------------
        function M = MapFcnMtx(this,r,s,t)
            % Isoparametric element shape
            M = this.ShpFcnMtx(r,s,t);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtx(this,r,s,t)
            GradNpar = this.gradShpFcnMtxVol(r,s,t);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxEdge(this,r)
             GradNpar = this.edge.gradShpFcnMtx(r,[],[]);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxFace(this,r,s)
            GradNpar = this.face.gradShpFcnMtx(r,s,[]);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxVol(this,r,s,t)
            GradNpar = zeros(this.dim,this.nen);
            
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