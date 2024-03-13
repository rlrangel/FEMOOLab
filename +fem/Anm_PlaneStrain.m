%% Anm_PlaneStrain Class (Plane Strain Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with plane strain models in a structural analysis.
%
%% Class definition
%
classdef Anm_PlaneStrain < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_PlaneStrain(meth)
            this = this@fem.Anm(fem.Anm.STRUCTURAL,fem.Anm.PLANE_STRAIN,meth,2,3,[1 2]);
            
            % Types of response
            this.DISPL_X  = true;  % Displacement X
            this.DISPL_Y  = true;  % Displacement Y
            this.SIGMA_XX = true;  % Normal stress XX
            this.SIGMA_YY = true;  % Normal stress YY
            this.SIGMA_ZZ = true;  % Normal stress ZZ
            this.TAU_XY   = true;  % Shear stress XY
            this.SIGMA_1  = true;  % Principal stress 1
            this.SIGMA_2  = true;  % Principal stress 2
            this.TAU_MAX  = true;  % Maximum shear stress
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anm
    methods
        %------------------------------------------------------------------
        % Set initial properties of elements.
        function setupElemProps(~,mdl)
            for i = 1:mdl.nel
                mdl.elems(i).TGNmtx(); % Gauss-to-node results transformation matrix
            end
        end
        
        %------------------------------------------------------------------
        % Assemble material constitutive matrix of a given element.
        function C = Cmtx(~,elem)
            E = elem.mat.E;
            v = elem.mat.v;
            e = E/((1+v)*(1-(2*v)));
            
            C = e * [ 1-v  v    0;
                      v    1-v  0;
                      0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble gradient matrix at a given position in parametric
        % coordinates of an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(this,elem,GradNcar,~,~)
            B = zeros(this.ndvc,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0;
                B(2,2*i-1) = 0;               B(2,2*i) = GradNcar(2,i);
                B(3,2*i-1) = GradNcar(2,i);   B(3,2*i) = GradNcar(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Return the ridigity coefficient at a given position in
        % parametric coordinates of an element.
        function coeff = rigidityCoeff(~,~,~,~)
            coeff = 1;
        end
        
        %------------------------------------------------------------------
        % Return the mass coefficient of an element.
        function coeff = massCoeff(~,elem)
            coeff = elem.mat.rho;
        end
        
        %------------------------------------------------------------------
        % Assemble global stiffness matrix.
        function K = gblStiffMtx(~,mdl)
            % Initialize global stiffness matrix
            K = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element matrices
                if mdl.anm.meth == mdl.anm.ISOPARAMETRIC
                    Kdiff = mdl.elems(i).stiffDiffMtx(); % diffusive term
                elseif mdl.anm.meth == mdl.anm.ISOGEOMETRIC
                    surface = mdl.surfaces(mdl.elems(i).surfId);
                    Kdiff = mdl.elems(i).stiffDiffMtx(surface);
                end
                
                % Assemble element matrices to global matrix
                gle = mdl.elems(i).gle;
                K(gle,gle) = K(gle,gle) + Kdiff;
            end
        end
        
        %------------------------------------------------------------------
        % Modify system arrays to include stabilization components for the
        % convective term.
        function [K,F] = stabConvec(~,~,K,F)
            return;
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to 1st time derivative of d.o.f.'s.
        % Damping matrix in structural analysis.
        function C = gblRate1Mtx(~,mdl)
            % NOT IMPLEMENTED
            C = zeros(mdl.neq,mdl.neq);
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to 2nd time derivative of d.o.f.'s.
        % Mass matrix in structural analysis.
        function M = gblRate2Mtx(~,mdl)
            % Initialize global mass matrix
            M = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element matrix
                Me = mdl.elems(i).massMtx();
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                M(gle,gle) = M(gle,gle) + Me;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        % Input:
        %  C: constituive matrix
        %  B: gradient matrix
        %  u: displacements results
        function str = pointDerivedVar(~,C,B,u)
            str = C * B * u;
        end
        
        %------------------------------------------------------------------
        % Compute derived variables and the principal values and directions
        % at Gauss points of all elements.
        % In structural analysis, derived variables are stresses.
        function gaussDerivedVar(this,mdl)
            r = mdl.res;
            
            npts = 0;
            for i = 1:mdl.nel
                % Compute stress components and cartesian coord at Gauss points
                U = r.U(mdl.elems(i).gle);
                
                if mdl.anm.meth == mdl.anm.ISOPARAMETRIC
                    [ngp,str,gpc] = mdl.elems(i).derivedVar(U);
                elseif mdl.anm.meth == mdl.anm.ISOGEOMETRIC
                    surface = mdl.surfaces(mdl.elems(i).surfId);
                    [ngp,str,gpc] = mdl.elems(i).derivedVar(U,surface);
                end
                
                r.ngp(i) = ngp;
                r.sxx_gp(1:ngp,i) = str(1,1:ngp);
                r.syy_gp(1:ngp,i) = str(2,1:ngp);
                r.txy_gp(1:ngp,i) = str(3,1:ngp);
                
                for j = 1:ngp
                    npts = npts + 1;
                    
                    % Store cartesian coordinates
                    r.x_gp(npts) = gpc(1,j);
                    r.y_gp(npts) = gpc(2,j);
                    
                    % Compute principal stresses
                    [prc,thetap]   = this.princStress(str(:,j));
                    r.s1_gp(j,i)   = prc(1);
                    r.s2_gp(j,i)   = prc(2);
                    r.tmax_gp(j,i) = prc(3);
                    
                    % Components of principal stresses in X-Y directions
                    r.s1x_gp(npts) = r.s1_gp(j,i)*cos(thetap);
                    r.s1y_gp(npts) = r.s1_gp(j,i)*sin(thetap);
                    r.s2x_gp(npts) = r.s2_gp(j,i)*cos(thetap+(pi/2.0));
                    r.s2y_gp(npts) = r.s2_gp(j,i)*sin(thetap+(pi/2.0));
                end
            end
        end
        
        %------------------------------------------------------------------
        % Extrapolate Gauss point results of derived variables to element
        % node results.
        % The nodal results are computed by extrapolation of Gauss point
        % results using the TGN matrix.
        % In structural analysis, derived variables are stresses.
        function elemDerivedVarExtrap(~,mdl)
            res = mdl.res;
            
            for i = 1:mdl.nel
                TGN = mdl.elems(i).TGN;
                nep = mdl.elems(i).shape.nep;
                ngp = res.ngp(i);
                
                gle = mdl.elems(i).gle;
                gle_x = gle(1:2:end);
                gle_y = gle(2:2:end);
                
                U_x = mdl.res.U(gle_x);
                U_y = mdl.res.U(gle_y);
                
                if mdl.anm.meth == mdl.anm.ISOPARAMETRIC
                    res.dx_elemextrap(:,i) = U_x;
                    res.dy_elemextrap(:,i) = U_y;
                    
                elseif mdl.anm.meth == mdl.anm.ISOGEOMETRIC
                    surfId = mdl.elems(i).surfId;
                    surface = mdl.surfaces(surfId);
                    xiSpan = mdl.elems(i).knotSpanXi;
                    etaSpan = mdl.elems(i).knotSpanEta;

                    % Surface properties
                    p = surface.degreeXi;
                    q = surface.degreeEta;
                    knotVectorXi = surface.knotVectorXi;
                    knotVectorEta = surface.knotVectorEta;
                    weights = surface.weights;

                    for j = 1:nep
                        % Extrapolation points parent coordinates
                        r = mdl.elems(i).shape.parCoord(j,1);
                        s = mdl.elems(i).shape.parCoord(j,2);

                        % Extrapolation points parametric coordinates
                        xi  = mdl.elems(i).shape.parent2ParametricSpace(xiSpan,r);
                        eta = mdl.elems(i).shape.parent2ParametricSpace(etaSpan,s);

                        % Basis functions
                        N = mdl.elems(i).shape.Nmtx(xi,eta,p,q,knotVectorXi,knotVectorEta,weights);

                        res.dx_elemextrap(j,i) = N * U_x;
                        res.dy_elemextrap(j,i) = N * U_y;
                    end
                end

                res.sxx_elemextrap(1:nep,i)  = TGN * res.sxx_gp(1:ngp,i);
                res.syy_elemextrap(1:nep,i)  = TGN * res.syy_gp(1:ngp,i);
                res.txy_elemextrap(1:nep,i)  = TGN * res.txy_gp(1:ngp,i);
                res.s1_elemextrap(1:nep,i)   = TGN * res.s1_gp(1:ngp,i);
                res.s2_elemextrap(1:nep,i)   = TGN * res.s2_gp(1:ngp,i);
                res.tmax_elemextrap(1:nep,i) = TGN * res.tmax_gp(1:ngp,i);
            end
        end
        
        %------------------------------------------------------------------
        % Smooth element node results of derived variables to global
        % node results.
        % The nodal global nodal results are computed by averaging values
        % of element extrapolated nodal results of all elements adjacent
        % to each node.
        % In structural analysis, derived variables are stresses.
        function nodeDerivedVarExtrap(~,mdl)
            res = mdl.res;
            
            % Sum contributions of extrapolated node results from connected elements
            for i = 1:mdl.nel
                for j = 1:mdl.elems(i).shape.nep
                    n = mdl.elems(i).shape.extNodes(j).id; % num da linha da matriz dos ePoints
                    res.dx_nodeextrap(n)   = res.dx_nodeextrap(n)   + res.dx_elemextrap(j,i);
                    res.dy_nodeextrap(n)   = res.dy_nodeextrap(n)   + res.dy_elemextrap(j,i);
                    res.sxx_nodeextrap(n)  = res.sxx_nodeextrap(n)  + res.sxx_elemextrap(j,i);
                    res.syy_nodeextrap(n)  = res.syy_nodeextrap(n)  + res.syy_elemextrap(j,i);
                    res.txy_nodeextrap(n)  = res.txy_nodeextrap(n)  + res.txy_elemextrap(j,i);
                    res.s1_nodeextrap(n)   = res.s1_nodeextrap(n)   + res.s1_elemextrap(j,i);
                    res.s2_nodeextrap(n)   = res.s2_nodeextrap(n)   + res.s2_elemextrap(j,i);
                    res.tmax_nodeextrap(n) = res.tmax_nodeextrap(n) + res.tmax_elemextrap(j,i);
                end
            end
            
            % Average nodal values by the number of connected elements
            for i = 1:mdl.nep
                nadjelems = length(mdl.extNodes(i).elems);
                res.dx_nodeextrap(i)   = res.dx_nodeextrap(i)   / nadjelems;
                res.dy_nodeextrap(i)   = res.dy_nodeextrap(i)   / nadjelems;
                res.sxx_nodeextrap(i)  = res.sxx_nodeextrap(i)  / nadjelems;
                res.syy_nodeextrap(i)  = res.syy_nodeextrap(i)  / nadjelems;
                res.txy_nodeextrap(i)  = res.txy_nodeextrap(i)  / nadjelems;
                res.s1_nodeextrap(i)   = res.s1_nodeextrap(i)   / nadjelems;
                res.s2_nodeextrap(i)   = res.s2_nodeextrap(i)   / nadjelems;
                res.tmax_nodeextrap(i) = res.tmax_nodeextrap(i) / nadjelems;
            end
        end

        %------------------------------------------------------------------
        function extPointsCoords(~,mdl)
            if mdl.anm.meth == mdl.anm.ISOPARAMETRIC
                for i = 1:mdl.nel
                    mdl.elems(i).shape.setExtNodesCoord();
                end
                mdl.nep = mdl.nnp;
                
                extNodes(mdl.nep,1) = fem.ExtNode();

                mdl.extNodes = extNodes;
                for i = 1:mdl.nep
                    mdl.extNodes(i).id = mdl.nodes(i).id;
                    mdl.extNodes(i).coord = mdl.nodes(i).coord;
                    mdl.extNodes(i).elems = fem.Element_Isoparametric().empty;
                    mdl.extNodes(i).elems = mdl.nodes(i).elems;
                end
                
                for i = 1:mdl.nel
                    for j = 1:mdl.elems(i).shape.nep
                        extNodeId = mdl.elems(i).shape.nodes(j).id;
                        mdl.elems(i).shape.extNodes(j) = mdl.extNodes(extNodeId);
                    end
                end
                
            elseif mdl.anm.meth == mdl.anm.ISOGEOMETRIC
                GlobalExtNodes = [];
                for i = 1:mdl.nel
                    % Element knot spans
                    xiSpan = mdl.elems(i).knotSpanXi;
                    etaSpan = mdl.elems(i).knotSpanEta;

                    surfId = mdl.elems(i).surfId;
                    surface = mdl.surfaces(surfId);
                    mdl.elems(i).shape.setExtNodesCoord(xiSpan,etaSpan,surface);

                    extCarCoord = mdl.elems(i).shape.extCarCoord;
                    GlobalExtNodes = uniquetol([GlobalExtNodes; extCarCoord],1e-5,'ByRows',true);
                end
                mdl.nep = size(GlobalExtNodes,1);
                
                %extNodes(mdl.nep,1) = fem.ExtNode(2*ones(mdl.nep,1));
                extNodes(mdl.nep,1) = fem.ExtNode();

                mdl.extNodes = extNodes;
                for i = 1:mdl.nep
                    mdl.extNodes(i).id = i;
                    mdl.extNodes(i).coord = GlobalExtNodes(i,:);
                    mdl.extNodes(i).coord(1,3) = 0;
                    mdl.extNodes(i).elems = fem.Element_Isogeometric().empty;
                end
                
                for i = 1:mdl.nel
                    extCarCoord = mdl.elems(i).shape.extCarCoord;
                    [~, extNodeIds] = ismembertol(extCarCoord,GlobalExtNodes,1e-5,'ByRows',true);
                    mdl.elems(i).shape.setccwExtNodeIds(extNodeIds);
                    
                    % Set extrapolation nodes incidence
                    for j = 1:mdl.elems(i).shape.nep
                        mdl.extNodes(extNodeIds(j)).elems(end+1) = mdl.elems(i);
                        mdl.elems(i).shape.extNodes(j) = mdl.extNodes(extNodeIds(j));
                    end
                end
            end
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Compute principal stress components and orientation for a given
        % stress tensor.
        % Input:
        %  str:  stress tensor (sx, sy, txy) stored in a column vector.
        % Output:
        %  prc:    principal stress components (s1, s2, taumax) stored in
        %          a column vector.
        %  thetap: angle of normal of principal stress plane w.r.t x axis
        %          (angle is returned in radians from 0 to 180 degrees).
        function [prc,thetap] = princStress(str)
            sx  = str(1);
            sy  = str(2);
            txy = str(3);
            
            center = (sx+sy)/2;
            deltas = (sx-sy)/2;
            radius = sqrt((deltas^2) + (txy*txy));
            
            prc = zeros(3,1);
            prc(1) = center + radius;   % s1
            prc(2) = center - radius;   % s2
            prc(3) = radius;            % taumax
            
            if (abs(deltas) > 0.0)
                thetap = 0.5 * atan2(txy,deltas);
            elseif (txy > 0.0)
                thetap = pi / 4.0;
            elseif (txy < 0.0)
                thetap = -pi / 4.0;
            else
                thetap = 0.0;
            end
            
            if( thetap < 0.0 )
                thetap = pi + thetap;
            end
        end
    end
end