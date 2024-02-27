%% Anm_ThickPlate Class (Thick Plate Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with thick (Mindlin) plate models in a structural analysis.
%
%% Class definition
%
classdef Anm_ThickPlate < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_ThickPlate()
            this = this@fem.Anm(fem.Anm.STRUCTURAL,fem.Anm.THICK_PLATE,3,5,[3 4 5]);
            
            % Types of response
            this.DISPL_Z     = true;  % Displacement Z
            this.ROTAT_X     = true;  % Rotation X
            this.ROTAT_Y     = true;  % Rotation Y
            this.SHEAR_XZ    = true;  % Shear force XZ
            this.SHEAR_YZ    = true;  % Shear force YZ
            this.MOMENT_XX   = true;  % Moment XX
            this.MOMENT_YY   = true;  % Moment YY
            this.MOMENT_XY   = true;  % Moment XY
            this.MOMENT_1    = true;  % Principal moment 1
            this.MOMENT_2    = true;  % Principal moment 2
            this.TORSION_MAX = true;  % Maximum torsion moment
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
            e = E/(1-(v^2));
            t = elem.thk;
           
            C = e * [ 1   v   0        0                0;
                      v   1   0        0                0;
                      0   0   (1-v)/2  0                0;
                      0   0   0        (6*(1-v))/(t^2)  0;
                      0   0   0        0                (6*(1-v))/(t^2)];
        end
        
        %------------------------------------------------------------------
        % Assemble gradient matrix at a given position in parametric
        % coordinates of an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(this,elem,GradNcar,r,s)                                   
            % Matrix of d.o.f.'s shape functions w.r.t. parametric coordinates
            N = elem.shape.Nmtx(r,s);
                        
            % Assemble gradient matrix
            B = zeros(this.ndvc,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,3*i-2) = 0;               B(1,3*i-1) = GradNcar(1,i);  B(1,3*i) = 0;
                B(2,3*i-2) = 0;               B(2,3*i-1) = 0;              B(2,3*i) = GradNcar(2,i);
                B(3,3*i-2) = 0;               B(3,3*i-1) = GradNcar(2,i);  B(3,3*i) = GradNcar(1,i);
                B(4,3*i-2) = -GradNcar(1,i);  B(4,3*i-1) = N(1,i);         B(4,3*i) = 0;
                B(5,3*i-2) = -GradNcar(2,i);  B(5,3*i-1) = 0;              B(5,3*i) = N(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Return the ridigity coefficient at a given position in
        % parametric coordinates of an element.
        function coeff = rigidityCoeff(~,elem,~,~)
            coeff = elem.thk^3/12;
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
        % Compute stress components (mxx, myy, mxy, qxz, qyz) at a given 
        % point of an element.
        % Input:
        %  C: constituive matrix
        %  B: gradient matrix
        %  u: displacements results
        function str = pointDerivedVar(~,C,B,u)
            str =  C * B * u;
        end
        
        %------------------------------------------------------------------
        % Compute derived variables and the principal values and directions
        % at Gauss points of all elements.
        % In plate analysis, derived variables are internal forces.
        function gaussDerivedVar(this,mdl)
            r = mdl.res;
            
            npts = 0;
            for i = 1:mdl.nel
                % Compute internal forces components and cartesian coord at Gauss points
                U = r.U(mdl.elems(i).gle);
                [ngp,str,gpc] = mdl.elems(i).derivedVar(U);
                
                r.ngp(i)          = ngp;
                r.mxx_gp(1:ngp,i) = str(1,1:ngp);
                r.myy_gp(1:ngp,i) = str(2,1:ngp);
                r.mxy_gp(1:ngp,i) = str(3,1:ngp);
                r.qxz_gp(1:ngp,i) = str(4,1:ngp);
                r.qyz_gp(1:ngp,i) = str(5,1:ngp);
                
                for j = 1:ngp
                    npts = npts + 1;
                    
                    % Store cartesian coordinates
                    r.x_gp(npts) = gpc(1,j);
                    r.y_gp(npts) = gpc(2,j);
                    
                    % Compute principal moments
                    [prc,thetap]     = this.princStress(str(:,j));
                    r.m1_gp(j,i)     = prc(1);
                    r.m2_gp(j,i)     = prc(2);
                    r.tormax_gp(j,i) = prc(3);
                    
                    % Components of principal stresses in X-Y directions
                    r.m1x_gp(npts) = r.m1_gp(j,i)*cos(thetap);
                    r.m1y_gp(npts) = r.m1_gp(j,i)*sin(thetap);
                    r.m2x_gp(npts) = r.m2_gp(j,i)*cos(thetap+(pi/2.0));
                    r.m2y_gp(npts) = r.m2_gp(j,i)*sin(thetap+(pi/2.0));
                end
            end
        end
        
        %------------------------------------------------------------------
        % Extrapolate Gauss point results of derived variables to element
        % node results.
        % The nodal results are computed by extrapolation of Gauss point
        % results using the TGN matrix.
        % In plate analysis, derived variables are internal forces.
        function elemDerivedVarExtrap(~,mdl)
            r = mdl.res;
            for i = 1:mdl.nel
                TGN = mdl.elems(i).TGN;
                nen = mdl.elems(i).shape.nen;
                ngp = r.ngp(i);
                
                r.mxx_elemextrap(1:nen,i)    = TGN * r.mxx_gp(1:ngp,i);
                r.myy_elemextrap(1:nen,i)    = TGN * r.myy_gp(1:ngp,i);
                r.mxy_elemextrap(1:nen,i)    = TGN * r.mxy_gp(1:ngp,i);
                r.m1_elemextrap(1:nen,i)     = TGN * r.m1_gp(1:ngp,i);
                r.m2_elemextrap(1:nen,i)     = TGN * r.m2_gp(1:ngp,i);
                r.tormax_elemextrap(1:nen,i) = TGN * r.tormax_gp(1:ngp,i);
                r.qxz_elemextrap(1:nen,i)    = TGN * r.qxz_gp(1:ngp,i);
                r.qyz_elemextrap(1:nen,i)    = TGN * r.qyz_gp(1:ngp,i);
            end
        end
        
        %------------------------------------------------------------------
        % Smooth element node results of derived variables to global
        % node results.
        % The nodal global nodal results are computed by averaging values
        % of element extrapolated nodal results of all elements adjacent
        % to each node.
        % In plate analysis, derived variables are internal forces.
        function nodeDerivedVarExtrap(~,mdl)
            r = mdl.res;
            
            % Sum contributions of extrapolated node results from connected elements
            for i = 1:mdl.nel
                for j = 1:mdl.elems(i).shape.nen
                    n = mdl.elems(i).shape.nodes(j).id;
                    r.mxx_nodeextrap(n)    = r.mxx_nodeextrap(n)    + r.mxx_elemextrap(j,i);
                    r.myy_nodeextrap(n)    = r.myy_nodeextrap(n)    + r.myy_elemextrap(j,i);
                    r.mxy_nodeextrap(n)    = r.mxy_nodeextrap(n)    + r.mxy_elemextrap(j,i);
                    r.m1_nodeextrap(n)     = r.m1_nodeextrap(n)     + r.m1_elemextrap(j,i);
                    r.m2_nodeextrap(n)     = r.m2_nodeextrap(n)     + r.m2_elemextrap(j,i);
                    r.tormax_nodeextrap(n) = r.tormax_nodeextrap(n) + r.tormax_elemextrap(j,i);
                    r.qxz_nodeextrap(n)    = r.qxz_nodeextrap(n)    + r.qxz_elemextrap(j,i);
                    r.qyz_nodeextrap(n)    = r.qyz_nodeextrap(n)    + r.qyz_elemextrap(j,i);
                end
            end
            
            % Average nodal values by the number of connected elements
            for i = 1:mdl.nnp
                nadjelems = length(mdl.nodes(i).elems);
                r.mxx_nodeextrap(i)    = r.mxx_nodeextrap(i)    / nadjelems;
                r.myy_nodeextrap(i)    = r.myy_nodeextrap(i)    / nadjelems;
                r.mxy_nodeextrap(i)    = r.mxy_nodeextrap(i)    / nadjelems;
                r.m1_nodeextrap(i)     = r.m1_nodeextrap(i)     / nadjelems;
                r.m2_nodeextrap(i)     = r.m2_nodeextrap(i)     / nadjelems;
                r.tormax_nodeextrap(i) = r.tormax_nodeextrap(i) / nadjelems;
                r.qxz_nodeextrap(i)    = r.qxz_nodeextrap(i)    / nadjelems;
                r.qyz_nodeextrap(i)    = r.qyz_nodeextrap(i)    / nadjelems;
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