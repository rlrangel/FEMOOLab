%% Shape_Quad8 Class
%
%% Description
%
% This is a sub-class of the <shape.html Shape> class for the
% implementation of the *Serendipity Quadratic Quadrilateral Shape*
% (8-noded isoparametric quadrilateral).
%
% <<shape_quadrilateral8.png>>
%
classdef Shape_Quad8 < Shape
    %% Constructor method
    methods
        function this = Shape_Quad8(nodes)
            this = this@Shape(Shape.QUAD8,2,8);
            
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
            this.nodeIdsLcl = [ 1; 5; 2; 6; 3; 7; 4; 8 ];
            
            % Vector of global node ids (ccw order)
            this.nodeIdsGbl = [ nodes(1).id;
                                nodes(5).id;
                                nodes(2).id;
                                nodes(6).id;
                                nodes(3).id;
                                nodes(7).id;
                                nodes(4).id;
                                nodes(8).id ];
            
            % Matrix of cartesian coordinates of nodes [X Y]
            this.carCoord = [ nodes(1).coord(1) nodes(1).coord(2);
                              nodes(2).coord(1) nodes(2).coord(2);
                              nodes(3).coord(1) nodes(3).coord(2);
                              nodes(4).coord(1) nodes(4).coord(2);
                              nodes(5).coord(1) nodes(5).coord(2);
                              nodes(6).coord(1) nodes(6).coord(2);
                              nodes(7).coord(1) nodes(7).coord(2);
                              nodes(8).coord(1) nodes(8).coord(2) ];
            
            % Matrix of parametric coordinates of nodes [r s]
            this.parCoord = [ -1.0 -1.0;
                               1.0 -1.0;
                               1.0  1.0;
                              -1.0  1.0;
                               0.0 -1.0;
                               1.0  0.0;
                               0.0  1.0;
                              -1.0  0.0 ];
            
            % Sub-components
            this.edge = Shape_Lin3();
            
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
            
            N(5) = 0.50*(1-r*r)*(1-s);
            N(6) = 0.50*(1+r)*(1-s*s);
            N(7) = 0.50*(1-r*r)*(1+s);
            N(8) = 0.50*(1-r)*(1-s*s);
            N(1) = 0.25*(1-r)*(1-s) - 0.50*N(8) - 0.50*N(5);
            N(2) = 0.25*(1+r)*(1-s) - 0.50*N(5) - 0.50*N(6);
            N(3) = 0.25*(1+r)*(1+s) - 0.50*N(6) - 0.50*N(7);
            N(4) = 0.25*(1-r)*(1+s) - 0.50*N(7) - 0.50*N(8);
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
            
            GradNpar(1,1) = (2*r - 2*r*s - s*s   + s) / 4;
            GradNpar(2,1) = (2*s - r*r   - 2*r*s + r) / 4;
            GradNpar(1,2) = (2*r - 2*r*s + s*s   - s) / 4;
            GradNpar(2,2) = (2*s - r*r   + 2*r*s - r) / 4;
            GradNpar(1,3) = (2*r + 2*r*s + s*s   + s) / 4;
            GradNpar(2,3) = (2*s + r*r   + 2*r*s + r) / 4;
            GradNpar(1,4) = (2*r + 2*r*s - s*s   - s) / 4;
            GradNpar(2,4) = (2*s + r*r   - 2*r*s - r) / 4;
            GradNpar(1,5) = r*s - r;
            GradNpar(2,5) = (r*r - 1) / 2;
            GradNpar(1,6) = (1 - s*s) / 2;
            GradNpar(2,6) = -s - r*s;
            GradNpar(1,7) = -r - r*s;
            GradNpar(2,7) = (1 - r*r) / 2;
            GradNpar(1,8) = (s*s - 1) / 2;
            GradNpar(2,8) = r*s - s;
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
        % n2. The local id of the mid node at the target edge, is returned
        % in the parameter mid.
        function [valid,n1,n2,mid] = edgeLocalIds(this,corner1,corner2)
            valid = false;
            n1    = 0;
            n2    = 0;
            mid   = 0;

            % Get ids of corner nodes and check for consistency
            for i = 1:4
                node = this.nodes(i).id;
                if node == corner1
                    n1 = i;
                    break;
                end
            end
            
            if n1 == 0
                return;
            end
            
            for i = 1:4
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
                    mid = 5;
                elseif n2 == 4
                    mid = 8;
                else
                    n1 = 0;
                    n2 = 0;
                    return; 
                end
            elseif n1 == 2
                if n2 == 1
                    mid = 5;
                elseif n2 == 3
                    mid = 6;
                else
                    n1 = 0;
                    n2 = 0;
                    return; 
                end
            elseif n1 == 3
                if n2 == 2
                    mid = 6;
                elseif n2 == 4
                    mid = 7;
                else
                    n1 = 0;
                    n2 = 0;
                    return; 
                end
            else % if n1 == 4
                if n2 == 1
                    mid = 8;
                elseif n2 == 3
                    mid = 7;
                else
                    n1 = 0;
                    n2 = 0;
                    return; 
                end
            end
            
            valid = true;
        end
    end
end