%% Shape_Tet10 Class
%
%% Description
%
% This is a sub-class of the <shape.html Shape> class for the
% implementation of the *Quadratic Thetahedral Shape*
% (10-noded isoparametric tetrahedron).
%
% <<shape_tetrahedron10.png>>
%
classdef Shape_Tet10 < Shape
    %% Constructor method
    methods
        function this = Shape_Tet10(nodes)
            this = this@Shape(Shape.TET10,3,10);
            
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
            this.nodeIdsLcl = [ 1; 5; 2; 6; 3; 7; 8; 9; 10; 4 ];
            
            % Vector of global node ids
            this.nodeIdsGbl = [ nodes(1).id;
                                nodes(5).id;
                                nodes(2).id;
                                nodes(6).id;
                                nodes(3).id;
                                nodes(7).id;
                                nodes(8).id;
                                nodes(9).id;
                                nodes(10).id;
                                nodes(4).id ];
            
            % Matrix of cartesian coordinates of nodes [X Y Z]
            this.carCoord = [ nodes(1).coord(1)  nodes(1).coord(2)  nodes(1).coord(3);
                              nodes(2).coord(1)  nodes(2).coord(2)  nodes(2).coord(3);
                              nodes(3).coord(1)  nodes(3).coord(2)  nodes(3).coord(3);
                              nodes(4).coord(1)  nodes(4).coord(2)  nodes(4).coord(3);
                              nodes(5).coord(1)  nodes(5).coord(2)  nodes(5).coord(3);
                              nodes(6).coord(1)  nodes(6).coord(2)  nodes(6).coord(3);
                              nodes(7).coord(1)  nodes(7).coord(2)  nodes(7).coord(3);
                              nodes(8).coord(1)  nodes(8).coord(2)  nodes(8).coord(3);
                              nodes(9).coord(1)  nodes(9).coord(2)  nodes(9).coord(3);
                              nodes(10).coord(1) nodes(10).coord(2) nodes(10).coord(3) ];
            
            % Matrix of parametric coordinates of nodes [r s t]
            this.parCoord = [ 0.0  0.0  0.0;
                              1.0  0.0  0.0;
                              0.0  1.0  0.0;
                              0.0  0.0  1.0;
                              0.5  0.0  0.0;
                              0.5  0.5  0.0;
                              0.0  0.5  0.0;
                              0.0  0.0  0.5;
                              0.5  0.0  0.5;
                              0.0  0.5  0.5 ];
            
            % Sub-components
            this.edge = Shape_Lin3();
            this.face = Shape_Tri6();
            
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