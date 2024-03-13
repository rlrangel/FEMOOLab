%% Shape_Tria6 Class (Quadratic Triangular Element Shape)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <shape.html Shape: element shape super-class> to deal
% with 6-noded isoparametric triangle (quadratic triangular) elements:
%                           s
%                           ^
%                           |
%                         3 +
%                           |\
%                           | \
%                         6 +  + 5
%                           |   \ TRIA6
%                           |    \
%                         1 +--+--+ 2  ----> r
%                              4
%
%% Class definition
%
classdef Shape_Tria6 < fem.Shape
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape_Tria6(nodes)
            this = this@fem.Shape(fem.Shape.TRIA6,2,6,6);
            this.nep = this.nen;
            
            if (nargin > 0)
                this.nodes = nodes;
                
                % Cartesian nodal coordiantes matrix [X Y]
                this.carCoord = [nodes(1).coord(1) nodes(1).coord(2);
                                 nodes(2).coord(1) nodes(2).coord(2);
                                 nodes(3).coord(1) nodes(3).coord(2);
                                 nodes(4).coord(1) nodes(4).coord(2);
                                 nodes(5).coord(1) nodes(5).coord(2);
                                 nodes(6).coord(1) nodes(6).coord(2)];
                
                % Parametric nodal coordinates matrix [r s]
                this.parCoord = [0.0  0.0;
                                 1.0  0.0;
                                 0.0  1.0;
                                 0.5  0.0;
                                 0.5  0.5;
                                 0.0  0.5];
                
                % Vector of local node ids in ccw order
                this.ccwLocalExtNodeIds = [1  4  2  5  3  6];
                
                % Vector of global node ids in ccw order
                this.ccwExtNodeIds = ...
                [ nodes(1).id  nodes(4).id  nodes(2).id ...
                  nodes(5).id  nodes(3).id  nodes(6).id];
                
                % Area
                x = this.carCoord(1:3,1);
                y = this.carCoord(1:3,2);
                this.size = polyarea(x,y);
                
                % Characteristic length
                this.Lchr = this.size^(1/2);
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Shape
    methods
        %------------------------------------------------------------------
        function setExtNodesCoord(this)
            this.extCarCoord = this.carCoord;
        end
        
        %------------------------------------------------------------------
        % Computes shape size (area).
        function s = size(this)
            x = this.carCoord(1:3,1);
            y = this.carCoord(1:3,2);
            s = polyarea(x,y);
        end
        
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
            
            N(1) = 1 - 3*r - 3*s + 4*r*s + 2*r*r + 2*s*s;
            N(2) = 2*r*r - r;
            N(3) = 2*s*s - s;
            N(4) = 4*r - 4*r*r - 4*r*s;
            N(5) = 4*r*s;
            N(6) = 4*s - 4*r*s - 4*s*s;
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge d.o.f. shape functions at a given
        % position in parametric coordinates.
        function N = NmtxEdge(~,r)
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
            
            GradNpar(1,1) =  4*r + 4*s - 3;
            GradNpar(2,1) =  4*r + 4*s - 3;
            GradNpar(1,2) =  4*r - 1;
            GradNpar(2,2) =  0;
            GradNpar(1,3) =  0;
            GradNpar(2,3) =  4*s - 1;
            GradNpar(1,4) =  4 - 8*r - 4*s;
            GradNpar(2,4) = -4*r;
            GradNpar(1,5) =  4*s;
            GradNpar(2,5) =  4*r;
            GradNpar(1,6) = -4*s;
            GradNpar(2,6) =  4 - 4*r - 8*s;
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradMpar = gradMmtxEdge(~,r)
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