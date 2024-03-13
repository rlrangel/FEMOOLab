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
classdef Shape_Isogeometric_Bezier_Extraction < fem.Shape
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Shape_Isogeometric_Bezier_Extraction(nodes)
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
        function setExtNodesCoord(this,Ce,surface)
            for i = 1:this.nep
                % Extrapolation points parent coordinates
                r = this.parCoord(i,1);
                s = this.parCoord(i,2);
                
                % Polynomial orders
                p = surface.degreeXi;
                q = surface.degreeEta;
                
                % Bernstein polynomials
                Bernstein = this.Nmtx(r,s,p,q);
                Bernstein = Bernstein';
                
                % Weight vector
                ctrlpts = this.nodes;
                weights = arrayfun(@(x) x.weight, ctrlpts);

                % Diagonal weight matrix
                W = diag(weights);

                % Weight function
                weightsb = Ce' * weights;
                Wfunc = weightsb' * Bernstein;
                
                % Basis functions
                R = 1/Wfunc * W * Ce * Bernstein;
                
                % Cartesian coordinates
                X = this.carCoord;
                
                % Extrapolation points cartesian coordinates
                this.extCarCoord(i,:) = R' * X;
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
        function M = Mmtx(this,xi,eta,p,q)
            M = this.Nmtx(xi,eta,p,q);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions at a given position in
        % parametric coordinates.
        function N = Nmtx(this,xi,eta,p,q)
            N = zeros(1, (p+1) * (q+1));
            
            for j=1:q+1
                for i=1:p+1
                    N((p+1)*(j-1)+i)= this.bernstein(p,i,xi) * this.bernstein(q,j,eta);
                end
            end
            
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge d.o.f. shape functions at a given
        % position in parametric coordinates.
        function N = NmtxEdge(this,xi,p)
            N = zeros((p+1),1);
            for i = 1:p+1
                N(i) = this.bernstein(p,i,xi);
            end
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        % Since this is an isoparametric element shape it returns the
        % evaluation of displacement shape functions derivatives. 
        function GradMpar = gradMmtx(this,xi,eta,p,q)
            GradMpar = this.gradNmtx(xi,eta,p,q);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of d.o.f. shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradNpar = gradNmtx(this,xi,eta,p,q)
            dNdxi = zeros(1, (p+1) * (q+1));
            dNdeta = zeros(1, (p+1) * (q+1));
            
            for j=1:q+1
                for i=1:p+1
                    dNdxi((p+1)*(j-1)+i) = 0.5 * p * (this.bernstein(p-1,i-1,xi) - this.bernstein(p- 1,i,xi)) * this.bernstein(q,j,eta);
                    dNdeta((p+1)*(j-1)+i) = this.bernstein(p,i,xi) * 0.5 * q * (this.bernstein(q-1,j- 1,eta) - this.bernstein(q-1,j,eta));
                end
            end
            
            GradNpar = [dNdxi; dNdeta];
        end
        
        function B = bernstein(this,p,a,xi)
            if p==0 && a==1
                B=1;
            elseif p==0 && a~=1
                B=0;
            else
                if a<1 || a>p+1
                    B=0;
                else
                    B1=this.bernstein(p-1,a,xi); 
                    B2=this.bernstein(p-1,a-1,xi); 
                    B=0.5*(1-xi)*B1+0.5*(1+xi)*B2;
                end
            end
            
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge geometry shape functions derivatives
        % w.r.t. parametric coordinates at a given position.
        function GradMpar = gradMmtxEdge(this,xi,p)
            dNdxi = zeros((p+1),1);
            for i = 1:p+1
                dNdxi(i) = 0.5 * p * (this.bernstein(p-1,i-1,xi) - this.bernstein(p- 1,i,xi));
            end
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