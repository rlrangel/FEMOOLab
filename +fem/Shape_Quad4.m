%% Shape_Quad4 Class (Bilinear Quadrilateral Element Shape)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <shape.html Shape: element shape super-class> to deal
% with 4-noded isoparametric quadrilateral (bilinear quadrilateral) elements:
%
%
%                         4 +---------------+ 3
%                           |       s       |
%                           |       ^       |
%                           |       |       |
%                           |       ----> r |
%                           |               |
%                           |     QUAD4     |
%                           |               |
%                         1 +---------------+ 2
%
%
%% Class definition
%
classdef Shape_Quad4 < fem.Shape
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape_Quad4(nodes)
            this = this@fem.Shape(fem.Shape.QUAD4,2,4,4);
            this.nep = this.nen;
            
            if (nargin > 0)
                this.nodes = nodes;
                
                % Cartesian nodal coordiantes matrix [X Y]
                this.carCoord = [nodes(1).coord(1) nodes(1).coord(2);
                                 nodes(2).coord(1) nodes(2).coord(2);
                                 nodes(3).coord(1) nodes(3).coord(2);
                                 nodes(4).coord(1) nodes(4).coord(2)];
                
                % Parametric nodal coordinates matrix [r s]
                this.parCoord = [-1 -1;
                                  1 -1;
                                  1  1;
                                 -1  1];
                
                % Vector of local node ids in ccw order
                this.ccwLocalExtNodeIds = [1  2  3  4];
                
                % Vector of global node ids in ccw order
                this.ccwExtNodeIds = ...
                [nodes(1).id  nodes(2).id  nodes(3).id  nodes(4).id];
                
                % Area
                x = this.carCoord(:,1);
                y = this.carCoord(:,2);
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
            
            N(1) = 0.25*(1-r)*(1-s);
            N(2) = 0.25*(1+r)*(1-s);
            N(3) = 0.25*(1+r)*(1+s);
            N(4) = 0.25*(1-r)*(1+s);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge d.o.f. shape functions at a given
        % position in parametric coordinates.
        function N = NmtxEdge(~,r)
            N = zeros(1,2);
            
            N(1) = 0.5*(1-r);
            N(2) = 0.5*(1+r);
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
        % Evaluate matrix of edge geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradMpar = gradMmtxEdge(~,~)
            GradMpar = zeros(1,2);
            
            GradMpar(1,1) = -0.5;
            GradMpar(1,2) =  0.5;
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