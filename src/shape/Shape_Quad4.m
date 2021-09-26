%% Shape_Quad4 Class
%
%% Description
%
% This is a sub-class of the <shape.html Shape> class for the
% implementation of the *Bilinear Quadrilateral Shape*
% (4-noded isoparametric quadrilateral).
%
% <<shape_quadrilateral4.png>>
%
classdef Shape_Quad4 < Shape
    %% Constructor method
    methods
        function this = Shape_Quad4(nodes)
            this = this@Shape(Shape.QUAD4,2,4);
            
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
            
            % Vector of local node ids (ccw order)
            this.nodeIdsLcl = [ 1; 2; 3; 4 ];
            
            % Vector of global node ids (ccw order)
            this.nodeIdsGbl = [ nodes(1).id;
                                nodes(2).id;
                                nodes(3).id;
                                nodes(4).id ];
            
            % Matrix of cartesian coordinates of nodes [X Y]
            this.carCoord = [ nodes(1).coord(1) nodes(1).coord(2);
                              nodes(2).coord(1) nodes(2).coord(2);
                              nodes(3).coord(1) nodes(3).coord(2);
                              nodes(4).coord(1) nodes(4).coord(2) ];
            
            % Matrix of parametric coordinates of nodes [r s]
            this.parCoord = [ -1.0 -1.0;
                               1.0 -1.0;
                               1.0  1.0;
                              -1.0  1.0 ];
            
            % Sub-components
            this.edge = Shape_Lin2();
            
            % Characteristic length
            x = this.carCoord(:,1);
            y = this.carCoord(:,2);
            area = polyarea(x,y);
            this.len = area^(1/2);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtx(this,r,s,~)
            N = this.ShpFcnMtxFace(r,s);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxEdge(this,r) 
            N = this.edge.ShpFcnMtx(r,[],[]);
        end
        
        %------------------------------------------------------------------
        function N = ShpFcnMtxFace(this,r,s)
            N = zeros(1,this.nen);
            
            N(1) = 0.25*(1-r)*(1-s);
            N(2) = 0.25*(1+r)*(1-s);
            N(3) = 0.25*(1+r)*(1+s);
            N(4) = 0.25*(1-r)*(1+s);
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
        function GradNpar = gradShpFcnMtx(this,r,s,~)
            GradNpar = this.gradShpFcnMtxFace(r,s);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxEdge(this,r)
             GradNpar = this.edge.gradShpFcnMtx(r,[],[]);
        end
        
        %------------------------------------------------------------------
        function GradNpar = gradShpFcnMtxFace(this,r,s)
            GradNpar = zeros(this.dim,this.nen);
            
            GradNpar(1,1) = -0.25*(1-s);
            GradNpar(2,1) = -0.25*(1-r);
            GradNpar(1,2) =  0.25*(1-s);
            GradNpar(2,2) = -0.25*(1+r);
            GradNpar(1,3) =  0.25*(1+s);
            GradNpar(2,3) =  0.25*(1+r);
            GradNpar(1,4) = -0.25*(1+s);
            GradNpar(2,4) =  0.25*(1-r);
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
        
        %------------------------------------------------------------------
        % Returns the local ids of the nodes of a shape edge for a given
        % pair of corner nodes global ids.
        % In case the given corner nodes global ids do not correspond to
        % a valid element shape edge, the returned value for the valid
        % parameter is false. Otherwise, the returned value is true.
        % The two corner nodes local ids are returned in parameters n1 and
        % n2. The returned mid value is zero.
        function [valid,n1,n2,mid] = edgeLocalIds(this,corner1,corner2)
            valid = false;
            n1    = 0;
            n2    = 0;
            mid   = 0;

            % Get ids of corner nodes and check for consistency
            for i = 1:this.nen
                node = this.nodes(i).id;
                if node == corner1
                    n1 = i;
                    break;
                end
            end
            
            if n1 == 0
                return;
            end
            
            for i = 1:this.nen
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
            
            % Check for corner node consistency
            if n1 == 1
                if n2 ~= 2 && n2 ~= 4
                    n1 = 0;
                    n2 = 0;
                    return;
                end
            elseif n1 == 4
                if n2 ~= 1 && n2 ~= 3
                    n1 = 0;
                    n2 = 0;
                    return;
                end
            else
                if n2 ~= n1+1 && n2 ~= n1-1
                    n1 = 0;
                    n2 = 0;
                    return;
                end
            end
            
            valid = true;
        end
    end
end