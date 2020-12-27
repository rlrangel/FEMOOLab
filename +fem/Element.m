%% Element Class
%
%% Description
%
% This is a class that specifies a generic finite element in the
% FEMOOLab program. This element may have any shape and analysis
% model, and may work with any Gauss quadrature.
% All the specific behaviors of an element are treated by other objects
% that are properties of an element:
% * anm: object of <anm.html Anm: analysis model class>
% * shape: object of <shape.html Shape: element shape class>
% * material: object of <material.html Material: material class>
% * gauss: object of <gauss.html Gauss: integration quadrature class>
%
%% Class definition
%
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        id    int32 = int32.empty;              % identification number
        type  int32 = int32.empty;              % element shape type
        gle   int32 = int32.empty;              % gather vector (stores element global d.o.f. numbers)
        anm   = [];                             % object of Anm class
        shape = [];                             % object of Shape class
        
        % Physical attributes
        mat fem.Material = fem.Material.empty;  % object of Material class
        thk double       = double.empty;        % thickness
        
        % Gauss integration quadrature
        gauss = [];                             % object of Gauss class
        gstiff_order  int32  = int32.empty;     % order of Gauss quadrature for stiffness matrix computation
        gstress_order int32  = int32.empty;     % order of Gauss quadrature for stress computation
        gstress_npts  int32  = int32.empty;     % number of gauss points for stress computation
        TGN           double = double.empty;    % transformation matrix of Gauss points results to nodal results
        
        % Loads
        lineLoad   double = double.empty;       % matrix of uniform line loads [corner1,corner2,qx,qy,qz]
        domainLoad double = double.empty;       % vector of uniform domain load components [px,py,pz]
        
        % Fluxes
        lineFlux    double = double.empty;      % matrix of uniform line heat fluxes [corner1,corner2,q]
        lineConvec  double = double.empty;      % matrix of uniform line heat convection [corner1,corner2,h,Tenv]
        domainFlux  double = double.empty;      % uniform domain heat flux [G]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Element()
            this.type = fem.Shape.GENERIC;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Set Gauss object and properties.
        function setGauss(this,gauss,gstiff_order,gstress_order)
            this.gauss = gauss;
            this.gstiff_order = gstiff_order;
            this.gstress_order = gstress_order;
            [this.gstress_npts,this.TGN] = this.TGNmtx();
        end
        
        %------------------------------------------------------------------
        % Compute transformation matrix of gauss-to-node results.
        % Refs.:
        % -Hinton & Campbell, "Local and Global Smoothing of Discontinous
        % Finite Element Functions using a Least Squares Method",
        % Int. J. Num. Meth. Engng., Vol. 8, pp. 461-480, 1974.
        % -Burnett, D.S., "Finite Element Analysis - From Concepts to Applications",
        % Addison-Wesley, 1987.
        % -Martha, L.F., "Notas de Aula do Curso CIV 2118 - Metodo dos Elementos
        % Finitos", 1994.
        function [ngp,TGN] = TGNmtx(this)
            % Get parametric coordinates of Gauss points
            [ngp,~,gp] = this.gauss.quadrature(this.gstress_order);
            
            if (this.gstress_order == 1)
                TGN = ones(this.shape.nen,1);
            else
                % Compute S matrix that defines the coefficients of the
                % smoothing stress plane that fits the Gauss point stress
                % values in a least square sence
                P = zeros(3,3);
                P(1,1) = ngp;
                P(1,2) = sum(gp(1,:));
                P(1,3) = sum(gp(2,:));
                P(2,1) = P(1,2);
                P(2,2) = gp(1,:) * gp(1,:)';
                P(2,3) = gp(1,:) * gp(2,:)';
                P(3,1) = P(1,3);
                P(3,2) = P(2,3);
                P(3,3) = gp(2,:) * gp(2,:)';
                
                Q = ones(3,ngp);
                Q(2,:) = gp(1,:);
                Q(3,:) = gp(2,:);
                
                S = P\Q;
                
                % Compute nodal stress evaluation matrix, which is obtained
                % using nodal parametric coordinates in the smoothing stress
                % plane equation
                E = ones(this.shape.nen,3);
                E(:,2) = this.shape.parCoord(:,1);
                E(:,3) = this.shape.parCoord(:,2);
                
                % Compute transformation matrix
                TGN = E * S;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stiffness matrix.
        function k = stiffMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element matrix
            k = zeros(nen*ndof,nen*ndof);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Material constitutive matrix
            C = this.anm.Cmtx(this);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gstiff_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Matrix of geometry shape functions derivatives
                % w.r.t. parametric coordinates
                GradMpar = this.shape.gradMmtx(r,s);
                
                % Jacobian matrix
                J = GradMpar * X;
                detJ = det(J);
                
                % Matrix of d.o.f. shape functions derivatives
                % w.r.t. parametric coordinates
                GradNpar = this.shape.gradNmtx(r,s);
                
                % Matrix of d.o.f. shape functions derivatives
                % w.r.t. cartesian coordinates
                GradNcar = J \ GradNpar;
                
                % Strain matrix 
                B = this.anm.Bmtx(this,GradNcar,r,s);
                
                % Get the rigidity coefficient at this integration point
                rigdtyCoeff = this.anm.rigidityCoeff(this,r,s);
                
                % Accumulate Gauss point contributions
                k = k + w(i) * rigdtyCoeff * detJ * B' * C * B;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stiffness matrix portion that accounts for convection B.C.
        function k = stiffConvecMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;
            
            % Initialize element matrix
            k = zeros(nen*ndof,nen*ndof);
            
            % Loop over element line convection edges
            for q = 1:size(this.lineConvec,1)
                % Global IDs of edge initial and final nodes
                corner1 = this.lineConvec(q,1);
                corner2 = this.lineConvec(q,2);
                
                % Get local IDs of edge nodes
                [valid,n1,n2,mid] = this.shape.edgeLocalIds(corner1,corner2);
                if ~valid
                    continue;
                end

                % Compute number of edge nodes and assemble vector of edge nodes ids
                if mid == 0
                    nedgen = 2;
                    edgLocIds = [n1,n2];
                else
                    nedgen = 3;
                    edgLocIds = [n1,n2,mid];
                end

                % Initialize edge matrix
                kline = zeros(nedgen*ndof,nedgen*ndof);
                
                % Edge nodes coordinates
                X = this.shape.carCoord(edgLocIds,:);
                
                % Convection coefficient
                h = this.lineConvec(q,3);
                
                % Gauss points and weights for integration on edge
                [ngp,w,gp] = this.gauss.lineQuadrature(this.gstiff_order);
                
                % Loop over edge Gauss integration points
                for i = 1:ngp
                    % Parametric coordinates
                    r = gp(1,i);
                    
                    % Edge d.o.f. shape functions matrix
                    N = this.shape.NmtxEdge(n1,n2,r);
                    
                    % Matrix of edge geometry shape functions derivatives
                    % w.r.t. parametric coordinates
                    GradMpar = this.shape.gradMmtxEdge(n1,n2,r);
                    
                    % Jacobian matrix
                    J = GradMpar * X;
                    detJ = sqrt(J(1)*J(1) + J(2)*J(2));
                    
                    % Accumulate Gauss point contributions
                    kline = kline + w(i) * detJ * h * (N' * N);
                end
                
                % Edge gather vector (stores local d.o.f.'s numbers)
                gledge = zeros(nedgen*ndof,1);
                m = 0;
                for i = 1:nedgen
                    for j = 1:ndof
                        m = m + 1;
                        gledge(m) = (edgLocIds(i)-1)*ndof + j;
                    end
                end
                
                % Assemble edge matrix to element matrix
                k(gledge,gledge) = k(gledge,gledge) + kline;
            end
        end
        
        %------------------------------------------------------------------
        % Compute equivalent nodal forcing vector for element line forces
        % (load or flux), distributed over edges of an element.
        % Loads are assumed in global direction (loc_gbl is not used)
        function f = edgeEquivForceVct(this,lineForce)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;
            
            % Initialize element forcing vector
            f = zeros(nen*ndof,1);
            
            % Loop over element line forces
            for q = 1:size(lineForce,1)
                % Global IDs of edge initial and final nodes
                corner1 = lineForce(q,1);
                corner2 = lineForce(q,2);
                
                % Get local IDs of edge nodes
                [valid,n1,n2,mid] = this.shape.edgeLocalIds(corner1,corner2);
                if ~valid
                    continue;
                end

                % Compute number of edge nodes and assemble vector of edge nodes ids
                if mid == 0
                    nedgen = 2;
                    edgLocIds = [n1,n2];
                else
                    nedgen = 3;
                    edgLocIds = [n1,n2,mid];
                end

                % Initialize edge forcing vector
                fline = zeros(nedgen*ndof,1);
                
                % Edge nodes coordinates
                X = this.shape.carCoord(edgLocIds,:);
                
                % Load components
                p = lineForce(q,3:end);
                
                % Gauss points and weights for integration on edge
                [ngp,w,gp] = this.gauss.lineQuadrature(this.gstiff_order);
                
                % Loop over edge Gauss integration points
                for i = 1:ngp
                    % Parametric coordinates
                    r = gp(1,i);
                    
                    % Edge d.o.f. shape functions matrix
                    N = this.shape.NmtxEdge(n1,n2,r);
                    
                    % Matrix of edge geometry shape functions derivatives
                    % w.r.t. parametric coordinates
                    GradMpar = this.shape.gradMmtxEdge(n1,n2,r);
                    
                    % Jacobian matrix
                    J = GradMpar * X;
                    detJ = sqrt(J(1)*J(1) + J(2)*J(2));
                    
                    % Accumulate Gauss point contributions
                    m = 0;
                    for j = 1:nedgen
                        for k = 1:ndof
                            m = m + 1;
                            fline(m) = fline(m) + w(i) * detJ * N(j) * p(k);
                        end
                    end
                end
                
                % Edge gather vector (stores local d.o.f.'s numbers)
                gledge = zeros(nedgen*ndof,1);
                m = 0;
                for i = 1:nedgen
                    for j = 1:ndof
                        m = m + 1;
                        gledge(m) = (edgLocIds(i)-1)*ndof + j;
                    end
                end
                
                % Assemble edge forcing vetor to element forcing vector
                f(gledge) = f(gledge) + fline;
            end
        end
        
        %------------------------------------------------------------------
        % Compute equivalent nodal forcing vector for element domain forces
        % (load or flux).
        function f = domainEquivForceVct(this,domainForce)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element forcing vector
            f = zeros(nen*ndof,1);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gstiff_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = this.shape.Nmtx(r,s);
                
                % Matrix of geometry shape functions derivatives w.r.t.
                % parametric coordinates
                GradMpar = this.shape.gradMmtx(r,s);
                
                % Jacobian matrix
                J = GradMpar * X;
                detJ = det(J);
                
                % Accumulate Gauss point contributions
                m = 0;
                for j = 1:nen
                    for k = 1:ndof
                        m = m + 1;
                        f(m) = f(m) + w(i) * detJ * N(j) * domainForce(k);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute stress components at gauss points and gauss point cartesian
        % coordinates for a given element.
        % Input:
        %  d: solution at element nodes
        % Output:
        %  ngp:  number of Gauss points for stress evaluation
        %  str:  stress components (sx,sy,txy) at each Gauss point
        %  gpc:  Gauss point cartesian coordinates array
        function [ngp,str,gpc] = gaussStress(this,d)
            % Get Gauss points and weights
            [ngp,~,gp] = this.gauss.quadrature(this.gstress_order);
            
            % Initialize stress component matrix and Gauss point coordinates
            str = zeros(3,ngp);
            if (this.anm.type == fem.Anm.PLANE_CONDUCTION) % CLEAN IT !!!!
                str = zeros(2,ngp);
            end
            
            gpc = zeros(2,ngp);
            
            % Cartesian coordinates
            X = this.shape.carCoord;
            
            % Material constituive matrix
            C = this.anm.Cmtx(this);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Geometry shape functions matrix evaluated at this point
                M = this.shape.Mmtx(r,s);
                
                % Matrix of geometry shape functions derivatives w.r.t.
                % parametric coordinates
                GradMpar = this.shape.gradMmtx(r,s);
                
                % Jacobian matrix
                J = GradMpar * X;
                
                % Matrix of displacement shape functions derivatives
                % w.r.t. parametric coordinates
                GradNpar = this.shape.gradNmtx(r,s);
                
                % Matrix of displacement shape functions derivatives
                % w.r.t. cartesian coordinates
                GradNcar = J \ GradNpar;
                
                % Strain-displacemen matrix 
                B = this.anm.Bmtx(this,GradNcar,r,s);
                
                % Gauss point stress components and cartesian coordinates
                str(:,i) = this.anm.pointStress(C,B,d);
                gpc(:,i) = M * X;
            end
        end
    end
end