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
            nno = length(nodes);
            this = this@fem.Shape(fem.Shape.ISOGe,2,nno);
            
            if (nargin > 0)
                this.nodes = nodes;
                this.nep = 8;
                
                % Cartesian control points coordiantes matrix [X Y] and
                % weights matrix
                X = zeros(nno,1);
                Y = zeros(nno,1);
                weights = zeros(nno,1);
                for i = 1:nno
                    X(i) = nodes(i).coord(1);
                    Y(i) = nodes(i).coord(2);
                    weights(i) = nodes(i).weight;
                end
                carCoord = [X, Y];
                this.carCoord = carCoord;
                this.weights = weights;
                
                % Parametric nodal coordinates matrix [r s]
                tol = 1e-10;
                this.parCoord = [ -1+tol -1+tol;
                                   1-tol -1+tol;
                                   1-tol  1-tol;
                                  -1+tol  1-tol;
                                   0     -1+tol;
                                   1-tol  0;
                                   0      1-tol;
                                  -1+tol  0 ];
                              
                % Set zeros at extrapolation points cartesian coordinates
                this.epCarCoord = zeros(this.nep,2);

                % Vector of local node ids in ccw order
                this.ccwLocalNodeIds = [ 1  5  2  6  3  7  4  8 ];
                
                % Vector of global node ids in ccw order
                this.ccwNodeIds = ...
                [ nodes(1).id  nodes(5).id  nodes(2).id  nodes(6).id ...
                  nodes(3).id  nodes(7).id  nodes(4).id  nodes(8).id];
                
                % Area
                x = this.carCoord(1:4,1);
                y = this.carCoord(1:4,2);
                this.size = polyarea(x,y);
                
                % Characteristic length
                this.Lchr = this.size^(1/2);
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Shape
    methods
        function setEPoints(this,xiSpan,etaSpan,surface)
            for i = 1:this.nep
                % Extrapolation points parent coordinates
                r = this.parCoord(i,1);
                s = this.parCoord(i,2);

                % Extrapolation points parametric coordinates
                xi  = this.parent2ParametricSpace(xiSpan,r);
                eta = this.parent2ParametricSpace(etaSpan,s);
                
                % Surface properties
                p = surface.degreeU;
                q = surface.degreeV;
                knotVectorU = surface.knotVectorU;
                knotVectorV = surface.knotVectorV;
                weights = surface.weights;
                
                % Basis functions
                N = this.Nmtx(xi,eta,p,q,knotVectorU,knotVectorV,weights);
                
                % Extrapolation points cartesian coordinates
                this.epCarCoord(i,:) = N * this.carCoord;
            end
        end
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions at a given position in
        % parametric coordinates.
        % Since this is an isoparametric element shape it returns the
        % evaluation of d.o.f. shape functions. 
        function M = Mmtx(this,xi,eta,p,q,knotVectorU,knotVectorV,weights)
            M = this.Nmtx(xi,eta,p,q,knotVectorU,knotVectorV,weights);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given position in
        % parametric coordinates.
        function N = Nmtx(~,xi,eta,p,q,knotVectorU,knotVectorV,weights)
            addpath ders/
            [N, ~, ~] = NURBS2DBasisDers([xi;eta],p,q,knotVectorU,knotVectorV,weights');
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
        function GradMpar = gradMmtx(this,xi,eta,p,q,knotVectorU,knotVectorV,weights)
            GradMpar = this.gradNmtx(xi,eta,p,q,knotVectorU,knotVectorV,weights);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradNpar = gradNmtx(~,xi,eta,p,q,knotVectorU,knotVectorV,weights)
            addpath ders/
            [~, dNdxi, dNdeta] = NURBS2DBasisDers([xi;eta],p,q,knotVectorU,knotVectorV,weights');
            GradNpar = [dNdxi; dNdeta];
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
    
        %------------------------------------------------------------------
        function xi = parent2ParametricSpace(~,range,pt)
            xi = 0.5 * ( (range(2) - range(1)) * pt + range(2) + range(1)); 
        end
    end
end