%% Anm_AxisymStress Class  (Axisymmetric Stress Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with axisymmetric stress models in a structural analysis.
%
%% Class definition
%
classdef Anm_AxisymStress < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_AxisymStress()
            this = this@fem.Anm(fem.Anm.STRUCTURAL,fem.Anm.AXISYM_STRESS,2,4,[1 2]);
            
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
            
            C = e * [ 1-v  v    v    0;
                      v    1-v  v    0;
                      v    v    1-v  0;
                      0    0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble gradient matrix at a given position in parametric
        % coordinates of an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(this,elem,GradNcar,r,s)                
            % Geometry and d.o.f. shape functions matrix evaluated at this point
            M = elem.shape.Mmtx(r,s);
            N = elem.shape.Nmtx(r,s);
            
            % Location of evaluation point (X coordinate is the radius in axisymm.)
            p = M * elem.shape.carCoord;
            radius = p(1);
            
            % Assemble gradient matrix
            B = zeros(this.ndvc,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0.0;
                B(2,2*i-1) = N(i)/radius;     B(2,2*i) = 0.0;
                B(3,2*i-1) = 0.0;             B(3,2*i) = GradNcar(2,i);
                B(4,2*i-1) = GradNcar(2,i);   B(4,2*i) = GradNcar(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Return the ridigity coefficient at a given position in
        % parametric coordinates of an element.
        function coeff = rigidityCoeff(~,elem,r,s)
            % Geometry shape functions matrix evaluated at this point
            M = elem.shape.Mmtx(r,s);
            
            % Location of evaluation point (X coordinate is the radius in axisymm.)
            % For axisymmetric analysis, the rigidity coefficient is the
            % radius at the given point.
            p = M * elem.shape.carCoord;
            coeff = p(1);
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
                Kdiff = mdl.elems(i).stiffDiffMtx(); % diffusive term
                
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
            % Compute point stress components
            str_raw = C * B * u;
            
            % Skip tangential stress component
            str(1) = str_raw(1);
            str(2) = str_raw(3);
            str(3) = str_raw(4);
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
                [ngp,str,gpc] = mdl.elems(i).derivedVar(U);
                
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
            r = mdl.res;
            
            for i = 1:mdl.nel
                TGN = mdl.elems(i).TGN;
                nen = mdl.elems(i).shape.nen;
                ngp = r.ngp(i);
                
                r.sxx_elemextrap(1:nen,i)  = TGN * r.sxx_gp(1:ngp,i);
                r.syy_elemextrap(1:nen,i)  = TGN * r.syy_gp(1:ngp,i);
                r.txy_elemextrap(1:nen,i)  = TGN * r.txy_gp(1:ngp,i);
                r.s1_elemextrap(1:nen,i)   = TGN * r.s1_gp(1:ngp,i);
                r.s2_elemextrap(1:nen,i)   = TGN * r.s2_gp(1:ngp,i);
                r.tmax_elemextrap(1:nen,i) = TGN * r.tmax_gp(1:ngp,i);
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
            r = mdl.res;
            
            % Sum contributions of extrapolated node results from connected elements
            for i = 1:mdl.nel
                for j = 1:mdl.elems(i).shape.nen
                    n = mdl.elems(i).shape.nodes(j).id;
                    r.sxx_nodeextrap(n)  = r.sxx_nodeextrap(n)  + r.sxx_elemextrap(j,i);
                    r.syy_nodeextrap(n)  = r.syy_nodeextrap(n)  + r.syy_elemextrap(j,i);
                    r.txy_nodeextrap(n)  = r.txy_nodeextrap(n)  + r.txy_elemextrap(j,i);
                    r.s1_nodeextrap(n)   = r.s1_nodeextrap(n)   + r.s1_elemextrap(j,i);
                    r.s2_nodeextrap(n)   = r.s2_nodeextrap(n)   + r.s2_elemextrap(j,i);
                    r.tmax_nodeextrap(n) = r.tmax_nodeextrap(n) + r.tmax_elemextrap(j,i);
                end
            end
            
            % Average nodal values by the number of connected elements
            for i = 1:mdl.nnp
                nadjelems = length(mdl.nodes(i).elems);
                r.sxx_nodeextrap(i)  = r.sxx_nodeextrap(i)  / nadjelems;
                r.syy_nodeextrap(i)  = r.syy_nodeextrap(i)  / nadjelems;
                r.txy_nodeextrap(i)  = r.txy_nodeextrap(i)  / nadjelems;
                r.s1_nodeextrap(i)   = r.s1_nodeextrap(i)   / nadjelems;
                r.s2_nodeextrap(i)   = r.s2_nodeextrap(i)   / nadjelems;
                r.tmax_nodeextrap(i) = r.tmax_nodeextrap(i) / nadjelems;
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