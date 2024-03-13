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
classdef Element_Isogeometric < fem.Element
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        surfId      int32  = int32.empty;       % surface identification number
        knotSpanXi  double = double.empty;      % knot span in xi direction
        knotSpanEta double = double.empty;      % knot span in eta direction
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
        % ### Same as Element_Isoparametric ###
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
        % ### Not implemented (same as Element_Isoparametric) ###
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
        %%
        %------------------------------------------------------------------
        % Compute gradient matrix in a position of element domain given by
        % parametric coordinates.
        function B = BmtxElem(this,J2,Rders,xi,eta)
            % Matrix of d.o.f. shape functions derivatives w.r.t.
            % cartesian coordinates
            GradNcar = J2 \ Rders;
            
            % Gradient matrix
            B = this.anm.Bmtx(this,GradNcar,xi,eta);
        end

        %------------------------------------------------------------------
        % Compute Jacobian matrix associated with the mapping from the 
        % parent domain to the parametric domain.
        function J1 = JmtxDomainParentToParametric(this)
            % Element knot spans
            xiSpan  = this.knotSpanXi;
            etaSpan = this.knotSpanEta;
            
            % Jacobian matrix
            J1 = [0.5*(xiSpan(2) - xiSpan(1)), 0.0; 
                 0.0, 0.5*(etaSpan(2) - etaSpan(1))];
        end

        %------------------------------------------------------------------
        % Compute Jacobian matrix associated with the mapping from the
        % parametric domain to the physical domain in a position given by
        % parametric coordinates.
        function [J2,Rders] = JmtxDomain(this,surface,xi,eta)
            % Surface properties
            p = surface.degreeXi;
            q = surface.degreeEta;
            knotVectorXi = surface.knotVectorXi;
            knotVectorEta = surface.knotVectorEta;
            weights = surface.weights;
            
            % Matrix of basis functions derivatives
            Rders = this.shape.gradMmtx(xi,eta,p,q,knotVectorXi,...
                knotVectorEta,weights);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord;
            
            % Jacobian matrix
            J2 = Rders * X;
        end

        %------------------------------------------------------------------
        % Compute Jacobian matrix associated with the 1D mapping from the 
        % parent domain to the parametric domain.
        function J1 = JmtxEdgeParentToParametric(~,Span)
            % Jacobian matrix
            J1 = [0.5 *(Span(2) - Span(1)), 0.0]; 
        end
        
        %------------------------------------------------------------------
        % Compute Jacobian matrix associated with the 1D mapping from the
        % parametric domain to the physical domain in a position given by a
        % parametric coordinate.
        function J2 = JmtxEdge(this,edgLocIds,p,knotVector,weights,xi)
            % Basis functions derivative
            Nder = this.shape.gradMmtxEdge(xi,p,knotVector,weights);
            
            % Cartesian coordinates matrix
            X = this.shape.carCoord(edgLocIds,:);
            
            % Jacobian matrix
            J2 = Nder * X;
        end
        
        %------------------------------------------------------------------
        % Compute number of edge nodes and assemble vector of edge nodes ids
        function [nedgen,edgLocIds,direc] = edgeIds(this,corner1,corner2,degreeXi,degreeEta)
            % Get local IDs of edge nodes
            [valid,edgLocIds,direc] = this.shape.edgeLocalIds(corner1,corner2,degreeXi,degreeEta);
            
            % Check for valid edge nodes
            if (~valid)
                edgLocIds = [];
                direc = "";
                return;
            end
            
            % Number of edge nodes and vector of edge nodes ids
            nedgen = length(edgLocIds);
        end
        
        %------------------------------------------------------------------
        % Assemble edge gather vector (stores local d.o.f.'s numbers).
        % ### Same as Element_Isoparametric ###
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
                % Parent coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Element knot spans
                xiSpan  = this.knotSpanXi;
                etaSpan = this.knotSpanEta;
                
                % Element parametric coordinates
                xi  = this.shape.parentToParametricSpace(xiSpan,r);
                eta = this.shape.parentToParametricSpace(etaSpan,s);
            
                % Jacobian matrix parent to parametric space
                J1 = this.JmtxDomainParentToParametric;
                
                % Jacobian matrix parametric to physical space
                [J2,Rders] = this.JmtxDomain(surface,xi,eta);
                
                % Gradient matrix 
                B = this.BmtxElem(J2,Rders,xi,eta);
                
                % Rigidity coefficient at integration point
                h = this.anm.rigidityCoeff(this,xi,eta);
                
                % Accumulate Gauss point contributions
                K = K + w(i) * B' * C * B * h * det(J1) * det(J2);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % convective term: [N]'[v]'[B]
        % ### Not implemented (same as Element_Isoparametric) ###
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
                J = this.JmtxDomain(r,s);
                
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
        % ### Not implemented (same as Element_Isoparametric) ###
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
                J = this.JmtxDomain(r,s);
                
                % Gradient matrix 
                B = this.BmtxElem(J,r,s);
                
                % Accumulate Gauss point contributions
                K = K + w(i) * b * B' * Vmtx * B * det(J);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % radiative boundary conditions over edges: [N]'[N]h
        % ### Not implemented (same as Element_Isoparametric) ###
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
        % ### Not implemented (same as Element_Isoparametric) ###
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
                J = this.JmtxDomain(r,s);
                                
                % Accumulate Gauss point contributions
                M = M + w(i) * (N' * N) * m * det(J);
            end
        end
        
        %------------------------------------------------------------------
        % Numerical integration of equivalent nodal forcing vector from
        % standard natural boundary conditions prescribed over edges: [N]'q
        function F = edgeEquivForceVct(this,mdl)
            ndof = this.anm.ndof;
            nen  = this.shape.nen;
            
            % Surface
            surface = mdl.surfaces(this.surfId);
            
            % Initialize element forcing vector
            F = zeros(nen*ndof,1);
            
            % Loop over edges with standard NBC
            for q = 1:size(this.lineNBC1,1)
                % Global IDs of edge initial and final nodes
                corner1 = this.lineNBC1(q,1);
                corner2 = this.lineNBC1(q,2);
                
                % Number of edge nodes and vector of edge nodes IDs
                [nedgen,edgLocIds,direc] = this.edgeIds(corner1,corner2,surface.degreeXi,surface.degreeEta);
                if (nedgen == 0)
                    continue;
                end
                
                % Initialize edge forcing vector
                Fedge = zeros(nedgen*ndof,1);
                
                % Forcing components
                p = this.lineNBC1(q,3:end);
                
                % Gauss points and weights for integration on edge
                [ngp,w,gp] = this.gauss.lineQuadrature(this.gsystem_order);
                
                % Knot span, degree and knot vector
                if direc == "xi"
                    Span = this.knotSpanXi;
                    degree = surface.degreeXi;
                    knotVector = surface.knotVectorXi;
                elseif direc == "eta"
                    Span = this.knotSpanEta;
                    degree = surface.degreeEta;
                    knotVector = surface.knotVectorEta;
                else
                    return
                end
                
                % Weights along the curve
                if direc == "xi"
                    if edgLocIds(1) == 1
                        ctrlptsIds = surface.ctrlNet(1,:);
                    else
                        ctrlptsIds = surface.ctrlNet(end,:);
                    end
                    
                elseif direc == "eta"
                    if edgLocIds(1) == 1
                        ctrlptsIds = surface.ctrlNet(:,1);
                    else
                        ctrlptsIds = surface.ctrlNet(:,end);
                    end
                    
                else
                    return
                end
                    
                ctrlpts = mdl.nodes(ctrlptsIds);
                weights = arrayfun(@(x) x.weight, ctrlpts);
            
                % Loop over edge Gauss integration points
                for i = 1:ngp
                    % Parent coordinates
                    r = gp(i);
                    
                    % Element parametric coordinates
                    xi = this.shape.parentToParametricSpace(Span,r);
            
                    % Edge d.o.f. shape functions matrix
                    N = this.shape.NmtxEdge(xi,degree,knotVector,weights);
                    
                    % Jacobian parent to parametric space
                    J1 = this.JmtxEdgeParentToParametric(Span);
                    detJ1 = J1(1);
                    
                    % Jacobian parametric to physical space
                    J2 = this.JmtxEdge(edgLocIds,degree,knotVector,weights,xi);
                    detJ2 = sqrt(J2(1)*J2(1) + J2(2)*J2(2));
                    
                    % Plate Hole Analytical solution
                    coordCtrlPts = this.shape.carCoord(edgLocIds,:);
                    X = N * coordCtrlPts;
                    
                    x = X(1,1);
                    y = X(1,2);

                    a = 1;
                    r = sqrt(x*x + y*y);
                    theta = atan(y/x);

                    c2t = cos(2*theta);
                    c4t = cos(4*theta);
                    s2t = sin(2*theta);
                    s4t = sin(4*theta);
                    fac1 = (a/r)^2;
                    fac2 = fac1*fac1;
                    
                    exact_stress(1) = (1-fac1*(1.5*c2t+c4t)+1.5*fac2*c4t);
                    exact_stress(2) = (-fac1*(0.5*c2t-c4t)-1.5*fac2*c4t);
                    exact_stress(3) = (-fac1*(0.5*s2t+s4t)+1.5*fac2*s4t);
                    
                    p_analytical = zeros(2,1);
                    if x >= -4.001 && x <= -3.999
                        p_analytical(1) = exact_stress(1) * p(1) * (-10);
                        p_analytical(2) = exact_stress(3) * p(2) * (-10);
                    elseif y >= 3.999 && y<= 4.001
                        p_analytical(1) = exact_stress(3) * p(1) * (10);
                        p_analytical(2) = exact_stress(2) * p(2) * (10);
                    end
                    
                    % Accumulate Gauss point contributions
                    m = 0;
                    for j = 1:nedgen
                        for k = 1:ndof
                            m = m + 1;
%                             Fedge(m) = Fedge(m) + w(i) * N(j) * p(k) * detJ1 * detJ2;
                            Fedge(m) = Fedge(m) + w(i) * N(j) * p_analytical(k) * detJ1 * detJ2;
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
        % ### Not implemented (same as Element_Isoparametric) ###
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
                J = this.JmtxDomain(r,s);
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
        % ### Not implemented (same as Element_Isoparametric) ###
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
                J = this.JmtxDomain(r,s);
                
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
            
            % Surface properties
            p = surface.degreeXi;
            q = surface.degreeEta;
            knotVectorXi = surface.knotVectorXi;
            knotVectorEta = surface.knotVectorEta;
            weights = surface.weights;
            
            % Element knot spans
            xiSpan  = this.knotSpanXi;
            etaSpan = this.knotSpanEta;
            
            % Material constituive matrix
            C = this.anm.Cmtx(this);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                
                % Parent coordinates
                r = gp(1,i);
                s = gp(2,i);

                % Element parametric coordinates
                xi  = this.shape.parentToParametricSpace(xiSpan,r);
                eta = this.shape.parentToParametricSpace(etaSpan,s);
                
                % Basis functions
                R = this.shape.Mmtx(xi,eta,p,q,knotVectorXi,knotVectorEta,weights);
                
                % Gradient matrix 
                [J,Rders] = this.JmtxDomain(surface,xi,eta);
                B = this.BmtxElem(J,Rders,xi,eta);
                
                % Gauss point stress components and cartesian coordinates
                dvar(:,i) = this.anm.pointDerivedVar(C,B,u);
                gpc(:,i) = R * X;
            end
        end
    end
end