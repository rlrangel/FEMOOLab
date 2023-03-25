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
classdef Element_Isogeometric < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        id        int32 = int32.empty;          % identification number
        surfId    int32 = int32.empty;          % surface identification number
        knotSpanU double = double.empty;        % knot span in u direction
        knotSpanV double = double.empty;        % knot span in v direction
        gle       int32 = int32.empty;          % gather vector (stores element global d.o.f. numbers)
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
        
        % Natural boundary conditions (uniformly distributed over edges)
        lineNBC1 double = double.empty;         % matrix of standard NBCs (constant values) [corner1,corner2,forcing_components]
        lineNBC2 double = double.empty;         % matrix of radiative NBCs (d.o.f. dependent values) [corner1,corner2,forcing_term]
        
        % Convection conditions
        avgV   double = double.empty            % cartesian components of average nodal velocities
        normV  double = double.empty            % norm of average nodal velocities vector
        peclet double = double.empty            % Peclet number
        alpha  double = double.empty            % stabilization coefficient alpha
        beta   double = double.empty            % stabilization coefficient beta
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Element_Isogeometric()
            return;
        end
    end
    
    %% Public methods
    methods
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
                this.TGN = ones(this.shape.nep,1);
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
                E = ones(this.shape.nep,3);
                E(:,2) = this.shape.parCoord(:,1);
                E(:,3) = this.shape.parCoord(:,2);
                
                % Compute transformation matrix
                this.TGN = E * S;
            end
        end
        
        %------------------------------------------------------------------
        % Initialize convection properties: velocities, Peclet number, and
        % stabilization coefficients.
        function convProps(this)
            nen = this.shape.nen;
            dim = this.shape.dim;
            
            % Asemble nodal velocity vector            
            V = zeros(nen,dim);
            for i = 1:nen
                if (~isempty(this.shape.nodes(i).convVel))
                    V(i,:) = this.shape.nodes(i).convVel(1:dim);
                end
            end
            
            % Average nodal velocity vector and its norm
            this.avgV = mean(V,1);
            this.normV = norm(this.avgV);
            
            % Peclet number
            this.peclet = (this.normV*this.shape.Lchr)/(2*this.mat.k);
            
            % Stabilization coefficients
            this.alpha = coth(this.peclet) - 1/this.peclet;
            this.beta  = (this.alpha*this.shape.Lchr)/(2*this.normV);
        end
        
        %------------------------------------------------------------------
        % Compute gradient matrix in a position of element domain given by
        % parametric coordinates.
        function B = BmtxElem(this,J2,surface,r,s)
            
            % Element knot spans
            xiSpan  = this.knotSpanU;
            etaSpan = this.knotSpanV;
            
            % Element parametric coordinates
            xi  = this.shape.parent2ParametricSpace(xiSpan,r);
            eta = this.shape.parent2ParametricSpace(etaSpan,s);
            
            % Surface properties
            p = surface.degreeU;
            q = surface.degreeV;
            knotVectorU = surface.knotVectorU;
            knotVectorV = surface.knotVectorV;
            weights = surface.weights;
            
            % Matrix of d.o.f. shape functions derivatives w.r.t.
            % parametric coordinates
            GradNpar = this.shape.gradNmtx(xi,eta,p,q,knotVectorU,knotVectorV,weights);
            
            % Matrix of d.o.f. shape functions derivatives w.r.t.
            % cartesian coordinates
            GradNcar = J2 \ GradNpar;
            
            % Gradient matrix
            B = this.anm.Bmtx(this,GradNcar,r,s);
        end

        %------------------------------------------------------------------
        % Compute Jacobian matrix in a position of element domain given by
        % parametric coordinates.
        function J = JmtxDomainPaToPa(this)
            % Element knot spans
            xiSpan  = this.knotSpanU;
            etaSpan = this.knotSpanV;
            
            % Jacobian matrix
            J = [0.5*(xiSpan(2) - xiSpan(1)), 0.0; 
                 0.0, 0.5*(etaSpan(2) - etaSpan(1))];
        end
        
        %------------------------------------------------------------------
        % Compute Jacobian matrix in a position of element domain given by
        % parametric coordinates.
        function J = JmtxDomainPatoPhy(this,surface,r,s)
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Element knot spans
            xiSpan  = this.knotSpanU;
            etaSpan = this.knotSpanV;
            
            % Element parametric coordinates
            xi  = this.shape.parent2ParametricSpace(xiSpan,r);
            eta = this.shape.parent2ParametricSpace(etaSpan,s);
            
            % Surface properties
            p = surface.degreeU;
            q = surface.degreeV;
            knotVectorU = surface.knotVectorU;
            knotVectorV = surface.knotVectorV;
            weights = surface.weights;
            
            % Matrix of geometry shape functions derivatives w.r.t.
            % parametric coordinates
            GradMpar = this.shape.gradMmtx(xi,eta,p,q,knotVectorU,knotVectorV,weights);
            
            % Jacobian matrix
            J = GradMpar * X;
        end
        
        %------------------------------------------------------------------
        % Compute Jacobian matrix in a position of element edge given by
        % parametric coordinates.
        function J = JmtxEdge(this,edgLocIds,r)
            % Cartesian coordinates matrix
            X = this.shape.carCoord(edgLocIds,:);
            
            % Matrix of edge geometry shape functions derivatives w.r.t.
            % parametric coordinates
            GradMpar = this.shape.gradMmtxEdge(r);
            
            % Jacobian matrix
            J = GradMpar * X;
        end
        
        %------------------------------------------------------------------
        % Compute number of edge nodes and assemble vector of edge nodes ids
        function [nedgen,edgLocIds] = edgeIds(this,corner1,corner2)
            % Get local IDs of edge nodes
            [valid,n1,n2,mid] = this.shape.edgeLocalIds(corner1,corner2);
            
            % Check for valid edge nodes
            if (~valid)
                nedgen = 0;
                edgLocIds = [];
                return;
            end
            
            % Number of edge nodes and vector of edge nodes ids
            if (mid == 0)
                nedgen = 2;
                edgLocIds = [n1,n2];
            else
                nedgen = 3;
                edgLocIds = [n1,n2,mid];
            end
        end
        
        %------------------------------------------------------------------
        % Assemble edge gather vector (stores local d.o.f.'s numbers).
        function gledge = gleEdgeVct(this,nedgen,edgLocIds)
            ndof = this.anm.ndof;
            gledge = zeros(nedgen*ndof,1);
            m = 0;
            for i = 1:nedgen
                for j = 1:ndof
                    m = m + 1;
                    gledge(m) = (edgLocIds(i)-1)*ndof + j;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % diffusive term: [B]'[C][B]h
        function K = stiffDiffMtx(this,surface)
            % Initialize element matrix
            ndof = this.anm.ndof;
            nen = this.shape.nen;
            K = zeros(nen*ndof,nen*ndof);
            
            % Material constitutive matrix
            C = this.anm.Cmtx(this);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Jacobian matrix parent to parametric space
                J1 = this.JmtxDomainPaToPa;
                
                % Jacobian matrix parametric to physical space
                J2 = this.JmtxDomainPatoPhy(surface,r,s);
                
                % Gradient matrix 
                B = this.BmtxElem(J2,surface,r,s);
                
                % Rigidity coefficient at integration point
                h = this.anm.rigidityCoeff(this,r,s);
                
                % Accumulate Gauss point contributions
                K = K + w(i) * B' * C * B * h * det(J1) * det(J2);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % convective term: [N]'[v]'[B]
        function K = stiffConvMtx(this)
            % Initialize element matrix
            ndof = this.anm.ndof;
            nen = this.shape.nen;
            K = zeros(nen*ndof,nen*ndof);
            
            % Average nodal velocity vector
            v = this.avgV;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = this.shape.Nmtx(r,s);
                
                % Jacobian matrix
                J = this.JmtxDomainPatoPhy(r,s);
                
                % Gradient matrix 
                B = this.BmtxElem(J,r,s);
                
                % Accumulate Gauss point contributions
                K = K + w(i) * N' * v * B * det(J);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % stabilization of convective term: b[B]'[V][B]
        % (steady-state analysis by SUPG method)
        function K = stiffStabMtx(this)
            % Initialize element matrix
            ndof = this.anm.ndof;
            nen = this.shape.nen;
            K = zeros(nen*ndof,nen*ndof);
            
            % Velocity matrix
            Vmtx = this.avgV' * this.avgV;
            
            % Stabilization coefficient
            b = this.beta;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Jacobian matrix
                J = this.JmtxDomainPatoPhy(r,s);
                
                % Gradient matrix 
                B = this.BmtxElem(J,r,s);
                
                % Accumulate Gauss point contributions
                K = K + w(i) * b * B' * Vmtx * B * det(J);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % radiative boundary conditions over edges: [N]'[N]h
        function K = stiffRadMtx(this)
            % Initialize element matrix
            ndof = this.anm.ndof;
            nen = this.shape.nen;
            K = zeros(nen*ndof,nen*ndof);
            
            % Loop over edges with radiative NBC
            for q = 1:size(this.lineNBC2,1)
                % Global IDs of edge initial and final nodes
                corner1 = this.lineNBC2(q,1);
                corner2 = this.lineNBC2(q,2);
                
                % Number of edge nodes and vector of edge nodes IDs
                [nedgen,edgLocIds] = this.edgeIds(corner1,corner2);
                if (nedgen == 0)
                    continue;
                end
                
                % Initialize edge matrix
                Kedge = zeros(nedgen*ndof,nedgen*ndof);
                
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
                    
                    % Jacobian matrix
                    J = this.JmtxEdge(edgLocIds,r);
                    detJ = sqrt(J(1)*J(1) + J(2)*J(2));
                    
                    % Accumulate Gauss point contributions
                    Kedge = Kedge + w(i) * (N' * N) * h * detJ;
                end
                
                % Edge gather vector
                gledge = this.gleEdgeVct(nedgen,edgLocIds);
                
                % Assemble edge matrix to element matrix
                K(gledge,gledge) = K(gledge,gledge) + Kedge;
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of mass matrix: [N]'[N]m
        function M = massMtx(this)
            % Initialize element matrix
            ndof = this.anm.ndof;
            nen = this.shape.nen;
            M = zeros(nen*ndof,nen*ndof);
            
            % Mass coefficient
            m = this.anm.massCoeff(this);
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = this.shape.Nmtx(r,s);
                
                % Jacobian matrix
                J = this.JmtxDomainPatoPhy(r,s);
                                
                % Accumulate Gauss point contributions
                M = M + w(i) * (N' * N) * m * det(J);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of equivalent nodal forcing vector from
        % standard natural boundary conditions prescribed over edges: [N]'q
        function F = edgeEquivForceVct(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;
            
            % Initialize element forcing vector
            F = zeros(nen*ndof,1);
            
            % Loop over edges with standard NBC
            for q = 1:size(this.lineNBC1,1)
                % Global IDs of edge initial and final nodes
                corner1 = this.lineNBC1(q,1);
                corner2 = this.lineNBC1(q,2);
                
                % Number of edge nodes and vector of edge nodes IDs
                [nedgen,edgLocIds] = this.edgeIds(corner1,corner2);
                if (nedgen == 0)
                    continue;
                end
                
                % Initialize edge forcing vector
                Fedge = zeros(nedgen*ndof,1);
                
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
                    
                    % Jacobian matrix
                    J = this.JmtxEdge(edgLocIds,r);
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
                
                % Edge gather vector
                gledge = this.gleEdgeVct(nedgen,edgLocIds);
                
                % Assemble edge forcing vetor to element forcing vector
                F(gledge) = F(gledge) + Fedge;
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of equivalent nodal forcing vector from
        % internal domain source: [N]'Q
        function F = domainEquivForceVct(this)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;

            % Initialize element forcing vector
            F = zeros(nen*ndof,1);
            
            % Internal forcing source
            Q = this.src;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = this.shape.Nmtx(r,s);
                
                % Jacobian matrix
                J = this.JmtxDomainPatoPhy(r,s);
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
        % Numerical integration of nodal forcing vector accounting for the
        % stabilization of convective term: b[B]'[v]'Q
        % (steady-state analysis by SUPG method)
        function F = domainStabForceVct(this)
            % Initialize element forcing vector
            ndof = this.anm.ndof;
            nen = this.shape.nen;
            F = zeros(nen*ndof,1);
            
            % Internal forcing source, avg velocity vector, stabilization coeff
            Q = this.src;
            v = this.avgV;
            b = this.beta;
            
            % Gauss points and weights
            [ngp,w,gp] = this.gauss.quadrature(this.gsystem_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Jacobian matrix
                J = this.JmtxDomainPatoPhy(r,s);
                
                % Gradient matrix 
                B = this.BmtxElem(J,r,s);
                
                % Accumulate Gauss point contributions
                F = F + w(i) * b * B' * v' * Q * det(J);
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
        function [ngp,dvar,gpc] = derivedVar(this,u,surface)
            % Get Gauss points and weights
            [ngp,~,gp] = this.gauss.quadrature(this.gderive_order);
            
            % Initialize derived variable components and Gauss point coordinates
            dvar = zeros(this.anm.ndvc,ngp); 
            gpc  = zeros(this.shape.dim,ngp);
            
            % Cartesian coordinates
            X = this.shape.carCoord;
            
            % Material constituive matrix
            C = this.anm.Cmtx(this);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                
                % Parent coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Element knot spans
                xiSpan  = this.knotSpanU;
                etaSpan = this.knotSpanV;

                % Element parametric coordinates
                xi  = this.shape.parent2ParametricSpace(xiSpan,r);
                eta = this.shape.parent2ParametricSpace(etaSpan,s);

                % Surface properties
                p = surface.degreeU;
                q = surface.degreeV;
                knotVectorU = surface.knotVectorU;
                knotVectorV = surface.knotVectorV;
                weights = surface.weights;
                
                % Geometry shape functions matrix evaluated at this point
                M = this.shape.Mmtx(xi,eta,p,q,knotVectorU,knotVectorV,weights);
                
                % Gradient matrix 
                J = this.JmtxDomainPatoPhy(surface,r,s);
                B = this.BmtxElem(J,surface,r,s);
                
                % Gauss point stress components and cartesian coordinates
                dvar(:,i) = this.anm.pointDerivedVar(C,B,u);
                gpc(:,i) = M * X;
            end
        end
    end
end