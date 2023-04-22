%% Shape_Isogeometric Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <shape.html Shape: element shape super-class> to deal
% with 8-noded isoparametric quadrilateral (Serendipity quadratic
% quadrilateral) elements:
%                                  CTRL_PT_8
%               CTRL_PT_7+             +             +CTRL_PT_9
%                               ---------------
%                              |               |
%                              |   CTRL_PT_5   |
%               CTRL_PT_4+     |       +       |     +CTRL_PT_6
%                              |               |
%                              |     ISOGe     |
%                               ---------------
%               CTRL_PT_1+             +             +CTRL_PT_3
%                                  CTRL_PT_2
%
%% Class definition
%
classdef Shape_Isogeometric < fem.Shape
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape_Isogeometric(nodes)
            this = this@fem.Shape(fem.Shape.ISOGe,2,length(nodes),8);
            
            if (nargin > 0)
                this.nodes = nodes;
                
                % Cartesian control points coordiantes matrix [X Y]
                for i = 1:this.nen
                    this.carCoord(i,:) = [nodes(i).coord(1) nodes(i).coord(2)];
                end
                
                % Parametric nodal coordinates matrix [r s]
                tol = 1e-8;
                this.parCoord = [-1+tol -1+tol;
                                  1-tol -1+tol;
                                  1-tol  1-tol;
                                 -1+tol  1-tol;
                                  0     -1+tol;
                                  1-tol  0;
                                  0      1-tol;
                                 -1+tol  0];

                % Vector of local extrapolation nodes ids in ccw order
                this.ccwLocalExtNodeIds = [1  5  2  6  3  7  4  8];
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Shape
    methods
        function setExtNodesCoord(this,xiSpan,etaSpan,surface)
            for i = 1:this.nep
                % Extrapolation points parent coordinates
                r = this.parCoord(i,1);
                s = this.parCoord(i,2);

                % Extrapolation points parametric coordinates
                xi  = this.parent2ParametricSpace(xiSpan,r);
                eta = this.parent2ParametricSpace(etaSpan,s);
                
                % Surface properties
                p = surface.degreeXi;
                q = surface.degreeEta;
                knotVectorXi = surface.knotVectorXi;
                knotVectorEta = surface.knotVectorEta;
                weights = surface.weights;
                
                % Basis functions
                N = this.Nmtx(xi,eta,p,q,knotVectorXi,knotVectorEta,weights);
                
                % Extrapolation points cartesian coordinates
                this.extCarCoord(i,:) = N * this.carCoord;
            end
        end
        
        %------------------------------------------------------------------
        function setccwExtNodeIds(this,extNodeIds)
            this.ccwExtNodeIds = ...
            [extNodeIds(1)  extNodeIds(5)  extNodeIds(2)  extNodeIds(6)...
             extNodeIds(3)  extNodeIds(7)  extNodeIds(4)  extNodeIds(8)];
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions at a given position in
        % parametric coordinates.
        % Since this is an isoparametric element shape it returns the
        % evaluation of d.o.f. shape functions. 
        function M = Mmtx(this,xi,eta,p,q,knotVectorXi,knotVectorEta,weights)
            M = this.Nmtx(xi,eta,p,q,knotVectorXi,knotVectorEta,weights);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given position in
        % parametric coordinates.
        function N = Nmtx(~,xi,eta,p,q,knotVectorXi,knotVectorEta,weights)
            addpath ders/ders2d
            [N, ~, ~] = NURBS2DBasisDers([xi;eta],p,q,knotVectorXi,knotVectorEta,weights');
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge d.o.f. shape functions at a given
        % position in parametric coordinates.
        function N = NmtxEdge(this,xi,p,knotVector)
            weights = zeros(this.nen,1);
            for i = 1:this.nen
                weights(i) = this.nodes(i).weight;
            end
            addpath ders/ders1d
            [N, ~] = NURBS1DBasisDers(xi,p,knotVector,weights);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        % Since this is an isoparametric element shape it returns the
        % evaluation of displacement shape functions derivatives. 
        function GradMpar = gradMmtx(this,xi,eta,p,q,knotVectorXi,knotVectorEta,weights)
            GradMpar = this.gradNmtx(xi,eta,p,q,knotVectorXi,knotVectorEta,weights);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradNpar = gradNmtx(~,xi,eta,p,q,knotVectorXi,knotVectorEta,weights)
            addpath ders/ders2d
            [~, dNdxi, dNdeta] = NURBS2DBasisDers([xi;eta],p,q,knotVectorXi,knotVectorEta,weights');
            GradNpar = [dNdxi; dNdeta];
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradMpar = gradMmtxEdge(this,xi,p,knotVector)
            weights = zeros(this.nen,1);
            for i = 1:this.nen
                weights(i) = this.nodes(i).weight;
            end
            addpath ders/ders1d
            [~, dNdxi] = NURBS1DBasisDers(xi,p,knotVector,weights);
            GradMpar = dNdxi;
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
        function [valid,edgLocIds,direc] = edgeLocalIds(this,corner1,corner2,degreeXi,degreeEta)
            valid = false;
            n1    = 0;
            n2    = 0;
            edgLocIds   = [];
            direc = "";

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
                return;
            end
           
            % Check for corner node consistency and get mid point id
            %nodesIds_Vector = 1:(degreeXi+1)*(degreeEta+1);
            %nodesIds_Matrix = reshape(nodesIds_Vector,degreeEta+1,degreeXi+1);
            %nodesIds_Matrix = flip(nodesIds_Matrix');
            
            nodesIds_Vector = 1:(degreeXi+1)*(degreeEta+1);
            nodesIds_Matrix = reshape(nodesIds_Vector,degreeXi+1,degreeEta+1)';
            nodesIds_Matrix = flipud(nodesIds_Matrix);
            
            [row1, column1] = find(nodesIds_Matrix==n1);
            [row2, column2] = find(nodesIds_Matrix==n2);
            
            if row1 == row2
                edgLocIds = nodesIds_Matrix(row1,:);
                direc = "xi";
            elseif column1 == column2 
                edgLocIds = nodesIds_Matrix(:,column1)';
                edgLocIds = flip(edgLocIds);
                direc = "eta";
            else
                edgLocIds = [];
                direc = "";
                return;
            end
            valid = true;
        end
    
        %------------------------------------------------------------------
        function xi = parent2ParametricSpace(~,range,pt)
            xi = 0.5 * ( (range(2) - range(1)) * pt + range(2) + range(1)); 
        end
    end
end