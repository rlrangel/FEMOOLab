%% Shape_Quad8 Class (Serendipity Quadratic Quadrilateral Element Shape)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <shape.html Shape: element shape super-class> to deal
% with 8-noded isoparametric quadrilateral (Serendipity quadratic
% quadrilateral) elements:
%                                   7
%                         4 +-------+-------+ 3
%                           |       s       |
%                           |       ^       |
%                           |       |       |
%                         8 +       ----> r + 6
%                           |               |
%                           |     QUAD8     |
%                           |               |
%                         1 +-------+-------+ 2
%                                   5
%
%% Class definition
%
classdef Shape_Quad8 < fem.Shape
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape_Quad8(nodes)
            this = this@fem.Shape(fem.Shape.QUAD8,8);
            
            if (nargin > 0)
                this.nodes = nodes;
                
                % Cartesian nodal coordiantes matrix [X Y]
                this.carCoord = [ nodes(1).coord(1) nodes(1).coord(2);
                                  nodes(2).coord(1) nodes(2).coord(2);
                                  nodes(3).coord(1) nodes(3).coord(2);
                                  nodes(4).coord(1) nodes(4).coord(2);
                                  nodes(5).coord(1) nodes(5).coord(2);
                                  nodes(6).coord(1) nodes(6).coord(2);
                                  nodes(7).coord(1) nodes(7).coord(2);
                                  nodes(8).coord(1) nodes(8).coord(2) ];
                
                % Parametric nodal coordinates matrix [r s]
                this.parCoord = [ -1 -1;
                                   1 -1;
                                   1  1;
                                  -1  1;
                                   0 -1;
                                   1  0;
                                   0  1;
                                  -1  0 ];
                
                % Vector of local node ids in ccw order
                this.ccwLocalNodeIds = [ 1  5  2  6  3  7  4  8 ];
                
                % Vector of global node ids in ccw order
                this.ccwNodeIds = ...
                [ nodes(1).id  nodes(5).id  nodes(2).id  nodes(6).id ...
                  nodes(3).id  nodes(7).id  nodes(4).id  nodes(8).id];
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Shape
    methods
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions at a given position in
        % parametric coordinates.
        % Since this is an isoparametric element shape it returns the
        % evaluation of d.o.f. shape functions. 
        function M = Mmtx(this,r,s)
            M = this.Nmtx(r,s);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given position in
        % parametric coordinates.
        function N = Nmtx(this,r,s)
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
        % Evaluate matrix of edge d.o.f. shape functions at a given
        % position in parametric coordinates.
        function N = NmtxEdge(~,~,~,r)
            N = zeros(1,3);
            
            N(3) = 1-r*r;
            N(1) = 0.5*(1-r) - 0.5*N(3);
            N(2) = 0.5*(1+r) - 0.5*N(3);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        % Since this is an isoparametric element shape it returns the
        % evaluation of displacement shape functions derivatives. 
        function GradMpar = gradMmtx(this,r,s)
            GradMpar = this.gradNmtx(r,s);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradNpar = gradNmtx(this,r,s)
            GradNpar = zeros(2,this.nen);
            
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
        % Evaluate matrix of edge geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        % The edge is defined by two corner nodes local ids (n1,n2).
        function GradMpar = gradMmtxEdge(~,~,~,r)
            GradMpar = zeros(1,3);
            
            GradMpar(1,1) = -0.5 + r;
            GradMpar(1,2) =  0.5 + r;
            GradMpar(1,3) = -2*r;
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