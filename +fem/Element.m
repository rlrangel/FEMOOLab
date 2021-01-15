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
% * gauss: object of <gauss.html Gauss: integration quadrature class>
% * material: object of <material.html Material: material class>
%
%% Class definition
%
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        id    int32 = int32.empty;              % identification number
        gle   int32 = int32.empty;              % gather vector (stores element global d.o.f. numbers)
        anm   = [];                             % object of Anm class
        shape = [];                             % object of Shape class
        gauss = [];                             % object of Gauss class
        
        % Gauss integration quadrature
        gsystem_order int32  = int32.empty;     % order of Gauss quadrature for computation of global system arrays
        gderive_order int32  = int32.empty;     % order of Gauss quadrature for computation of derived variables
        gderive_npts  int32  = int32.empty;     % number of Gauss points for computation of derived variables        
        TGN           double = double.empty;    % transformation matrix of Gauss points results to nodal results
        
        % Physical attributes
        mat fem.Material = fem.Material.empty;  % object of Material class
        thk double       = double.empty;        % thickness
        
        % Domain forcing source
        src double = double.empty;              % components of body force (structural) / internal heat generation (thermal)
        
        % Natural boundary conditions
        lineNBC1 double = double.empty;         % matrix of uniform distributed values of natural boundary conditions over edges - constant values
                                                % [corner1,corner2,forcing_components]
        lineNBC2 double = double.empty;         % matrix of uniform distributed values of natural boundary conditions over edges - d.o.f. dependent values
                                                % [corner1,corner2,forcing_term]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Element()
            return;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute stiffness matrix part accounting for diffusion.
        function K = stiffDiffMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element matrix
            K = zeros(nen*ndof,nen*ndof);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Material constitutive matrix
            C = this.anm.Cmtx(this);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
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
                
                % Rigidity coefficient at integration point
                h = this.anm.rigidityCoeff(this,r,s);
                
                % Accumulate Gauss point contributions
                K = K + w(i) * B' * C * B * h * detJ;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stiffness matrix part accounting for convection.
        function K = stiffConvMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element matrix
            K = zeros(nen*ndof,nen*ndof);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Asemble nodal velocity vector (currently, X and Y componenets only!)
            V = zeros(nen,2);
            for i = 1:nen
                if (~isempty(this.shape.nodes(i).convVel))
                    V(i,1:2) = this.shape.nodes(i).convVel(1:2);
                end
            end
            
            % Average nodal velocity assumed for the entire element
            v = mean(V,1);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = this.shape.Nmtx(r,s);
                
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
                
                % Velocity at Gauss point
                %v = N * V;
                
                % Accumulate Gauss point contributions
                K = K + w(i) * N' * v * B * detJ;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stiffness matrix part accounting for radiation boundary
        % conditions (d.o.f. dependent natural B.C.'s) over edges.
        function K = stiffRadMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;
            
            % Initialize element matrix
            K = zeros(nen*ndof,nen*ndof);
            
            % Loop over edges with d.o.f. dependent NBC
            for q = 1:size(this.lineNBC2,1)
                % Global IDs of edge initial and final nodes
                corner1 = this.lineNBC2(q,1);
                corner2 = this.lineNBC2(q,2);
                
                % Get local IDs of edge nodes
                [valid,n1,n2,mid] = this.shape.edgeLocalIds(corner1,corner2);
                if (~valid)
                    continue;
                end
                
                % Compute number of edge nodes and assemble vector of edge nodes ids
                if (mid == 0)
                    nedgen = 2;
                    edgLocIds = [n1,n2];
                else
                    nedgen = 3;
                    edgLocIds = [n1,n2,mid];
                end
                
                % Initialize edge matrix
                Kedge = zeros(nedgen*ndof,nedgen*ndof);
                
                % Edge nodes coordinates
                X = this.shape.carCoord(edgLocIds,:);
                
                % Radiation coefficient
                h = this.lineNBC2(q,3);
                
                % Gauss points and weights for integration on edge
                [ngp,w,gp] = this.gauss.lineQuadrature(this.gsystem_order);
                
                % Loop over edge Gauss integration points
                for i = 1:ngp
                    % Parametric coordinates
                    r = gp(i);
                    
                    % Edge d.o.f. shape functions matrix
                    N = this.shape.NmtxEdge(r);
                    
                    % Matrix of edge geometry shape functions derivatives
                    % w.r.t. parametric coordinates
                    GradMpar = this.shape.gradMmtxEdge(r);
                    
                    % Jacobian matrix
                    J = GradMpar * X;
                    detJ = sqrt(J(1)*J(1) + J(2)*J(2));
                    
                    % Accumulate Gauss point contributions
                    Kedge = Kedge + w(i) * (N' * N) * h * detJ;
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
                K(gledge,gledge) = K(gledge,gledge) + Kedge;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stabilization matrix of convective term of steady-state
        % analysis by SUPG method.
        function K = stiffStabMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element matrix
            K = zeros(nen*ndof,nen*ndof);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Asemble nodal velocity vector (currently, X and Y componenets only!)
            V = zeros(nen,2);
            for i = 1:nen
                if (~isempty(this.shape.nodes(i).convVel))
                    V(i,1:2) = this.shape.nodes(i).convVel(1:2);
                end
            end
            
            % Average nodal velocity assumed for the entire element
            v = mean(V,1);
            
            % Velocity norm and matrix
            normV = norm(v);
            Vmtx = v' * v;
            
            % Characteristic length (assuming 2D)
            dim = 2;
            L = this.shape.size^(1/dim);
            
            % Peclet number
            Pe = (normV*L)/(2*this.mat.k);
            
            % Stabilization coefficients
            alpha = coth(abs(Pe)) - 1/abs(Pe);
            beta  = (alpha*L)/(2*normV);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
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
                
                % Accumulate Gauss point contributions
                K = K + w(i) * beta * B' * Vmtx * B * detJ;
            end
        end
        
        %------------------------------------------------------------------
        % Compute mass/capacity matrix.
        function M = massMtx(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element matrix
            M = zeros(nen*ndof,nen*ndof);
            
            % Mass coefficient
            m = this.anm.massCoeff(this);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = this.shape.Nmtx(r,s);
                
                % Matrix of geometry shape functions derivatives
                % w.r.t. parametric coordinates
                GradMpar = this.shape.gradMmtx(r,s);
                
                % Jacobian matrix
                J = GradMpar * X;
                detJ = det(J);
                                
                % Accumulate Gauss point contributions
                M = M + w(i) * (N' * N) * m * detJ;
            end
        end
        
        %------------------------------------------------------------------
        % Compute equivalent nodal forcing vector for constant natural
        % boundary conditions prescribed over edges.
        function F = edgeEquivForceVct(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;
            
            % Initialize element forcing vector
            F = zeros(nen*ndof,1);
            
            % Loop over edges with constant NBC
            for q = 1:size(this.lineNBC1,1)
                % Global IDs of edge initial and final nodes
                corner1 = this.lineNBC1(q,1);
                corner2 = this.lineNBC1(q,2);
                
                % Get local IDs of edge nodes
                [valid,n1,n2,mid] = this.shape.edgeLocalIds(corner1,corner2);
                if (~valid)
                    continue;
                end

                % Compute number of edge nodes and assemble vector of edge nodes ids
                if (mid == 0)
                    nedgen = 2;
                    edgLocIds = [n1,n2];
                else
                    nedgen = 3;
                    edgLocIds = [n1,n2,mid];
                end

                % Initialize edge forcing vector
                Fedge = zeros(nedgen*ndof,1);
                
                % Edge nodes coordinates
                X = this.shape.carCoord(edgLocIds,:);
                
                % Forcing components
                p = this.lineNBC1(q,3:end);
                
                % Gauss points and weights for integration on edge
                [ngp,w,gp] = this.gauss.lineQuadrature(this.gsystem_order);
                
                % Loop over edge Gauss integration points
                for i = 1:ngp
                    % Parametric coordinates
                    r = gp(i);
                    
                    % Edge d.o.f. shape functions matrix
                    N = this.shape.NmtxEdge(r);
                    
                    % Matrix of edge geometry shape functions derivatives
                    % w.r.t. parametric coordinates
                    GradMpar = this.shape.gradMmtxEdge(r);
                    
                    % Jacobian matrix
                    J = GradMpar * X;
                    detJ = sqrt(J(1)*J(1) + J(2)*J(2));
                    
                    % Accumulate Gauss point contributions
                    m = 0;
                    for j = 1:nedgen
                        for k = 1:ndof
                            m = m + 1;
                            Fedge(m) = Fedge(m) + w(i) * N(j) * p(k) * detJ;
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
                F(gledge) = F(gledge) + Fedge;
            end
        end
        
        %------------------------------------------------------------------
        % Compute equivalent nodal forcing vector for internal domain source.
        function F = domainEquivForceVct(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element forcing vector
            F = zeros(nen*ndof,1);
            
            % Internal forcing source
            Q = this.src;
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
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
                        F(m) = F(m) + w(i) * N(j) * Q(k) * detJ;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute stabilization vector for internal domain source.
        function F = domainStabForceVct(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element forcing vector
            F = zeros(nen*ndof,1);
            
            % Internal forcing source
            Q = this.src;
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Asemble nodal velocity vector (currently, X and Y componenets only!)
            V = zeros(nen,2);
            for i = 1:nen
                if (~isempty(this.shape.nodes(i).convVel))
                    V(i,1:2) = this.shape.nodes(i).convVel(1:2);
                end
            end
            
            % Average nodal velocity assumed for the entire element
            v = mean(V,1);
            
            % Velocity norm
            normV = norm(v);
            
            % Characteristic length (assuming 2D)
            dim = 2;
            L = this.shape.size^(1/dim);
            
            % Peclet number
            Pe = (normV*L)/(2*this.mat.k);
            
            % Stabilization coefficients
            alpha = coth(abs(Pe)) - 1/abs(Pe);
            beta  = (alpha*L)/(2*normV);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
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
                
                % Accumulate Gauss point contributions
                F = F + w(i) * beta * B' * v' * Q * detJ;
            end
        end
        
        %------------------------------------------------------------------
        % Compute derived variables at Gauss points and gauss point cartesian
        % coordinates for a given element.
        % Input:
        %  u: solution at element nodes
        % Output:
        %  ngp:  number of Gauss points for derived variables evaluation
        %  dvar: derived variable components at each Gauss point
        %  gpc:  Gauss point cartesian coordinates array
        function [ngp,dvar,gpc] = derivedVar(this,u)
            % Get Gauss points and weights
            [ngp,~,gp] = this.gauss.quadrature(this.gderive_order);
            
            % Initialize derived variable components and Gauss point coordinates
            dvar = zeros(this.anm.ndvc,ngp); 
            gpc  = zeros(2,ngp);
            
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
                dvar(:,i) = this.anm.pointDerivedVar(C,B,u);
                gpc(:,i) = M * X;
            end
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
        function TGNmtx(this)
            % Get parametric coordinates of Gauss points
            [ngp,~,gp] = this.gauss.quadrature(this.gderive_order);
            
            % Set number of Gauss points for derived variables computation
            this.gderive_npts = ngp;
            
            if (this.gderive_order == 1)
                this.TGN = ones(this.shape.nen,1);
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
                this.TGN = E * S;
            end
        end
    end
end