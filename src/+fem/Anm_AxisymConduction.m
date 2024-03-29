%% Anm_AxisymConduction Class  (Axisymmetric Heat Conduction Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with axisymmetric heat conduction models in a thermal analysis.
%
%% Class definition
%
classdef Anm_AxisymConduction < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_AxisymConduction()
            this = this@fem.Anm(fem.Anm.THERMAL,fem.Anm.AXISYM_CONDUCTION,1,2,1);
            
            % Types of response
            this.TEMPERATURE = true;  % Temperature
            this.FLUX_XX     = true;  % flux XX
            this.FLUX_YY     = true;  % flux YY
            this.FLUX_MOD    = true;  % flux module
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anm
    methods
        %------------------------------------------------------------------
        % Set initial properties of elements.
        function setupElemProps(~,mdl)
            for i = 1:mdl.nel
                mdl.elems(i).TGNmtx();    % Gauss-to-node results transformation matrix
                mdl.elems(i).convProps(); % Convection properties
            end
        end
        
        %------------------------------------------------------------------
        % Assemble material constitutive matrix of a given element.
        function C = Cmtx(~,elem)
            k = elem.mat.k;
            
            C = [ k   0;
                  0   k ];
        end
        
        %------------------------------------------------------------------
        % Assemble gradient matrix at a given position in parametric
        % coordinates of an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(this,elem,GradNcar,~,~)                
            B = zeros(this.ndvc,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,i) = GradNcar(1,i);
                B(2,i) = GradNcar(2,i);
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
            coeff = elem.mat.rho * elem.mat.cp;
        end
        
        %------------------------------------------------------------------
        % Assemble global stiffness matrix.
        function K = gblStiffMtx(~,mdl)
            % Initialize global stiffness matrix
            K = zeros(mdl.neq,mdl.neq);
            
            % Get element matrices and assemble global matrix
            for i = 1:mdl.nel
                % Get element matrices
                Kdiff = mdl.elems(i).stiffDiffMtx(); % diffusive term
                Krad  = mdl.elems(i).stiffRadMtx();  % radiative B.C.
                
                % Assemble element matrices to global matrix
                gle = mdl.elems(i).gle;
                K(gle,gle) = K(gle,gle) + Kdiff + Krad;
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
        % Capacity matrix in thermal analysis.
        function C = gblRate1Mtx(~,mdl)
            % Initialize global capacity matrix
            C = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element matrix
                Ce = mdl.elems(i).massMtx();
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                C(gle,gle) = C(gle,gle) + Ce;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to 2nd time derivative of d.o.f.'s.
        function M = gblRate2Mtx(~,mdl)
            M = zeros(mdl.neq,mdl.neq);
        end
        
        %------------------------------------------------------------------
        % Compute flux components (fx, fy) at a given point of an element.
        % Input:
        %  C: constituive matrix
        %  B: gradient matrix
        %  u: temperature results
        function flx = pointDerivedVar(~,C,B,u)
            flx = -C * B * u;
        end
        
        %------------------------------------------------------------------
        % Compute derived variables and the principal values and directions
        % at Gauss points of all elements.
        % In thermal analysis, derived variables are heat fluxes.
        function gaussDerivedVar(~,mdl)
            r = mdl.res;
            
            npts = 0;
            for i = 1:mdl.nel
                % Compute flux components and cartesian coord at Gauss points
                U = r.U(mdl.elems(i).gle);
                [ngp,flx,gpc] = mdl.elems(i).derivedVar(U);
                
                r.ngp(i) = ngp;
                r.fxx_gp(1:ngp,i) = flx(1,1:ngp);
                r.fyy_gp(1:ngp,i) = flx(2,1:ngp);
                
                for j = 1:ngp
                    npts = npts + 1;
                    
                    % Store cartesian coordinates
                    r.x_gp(npts) = gpc(1,j);
                    r.y_gp(npts) = gpc(2,j);
                    
                    % Compute flux module
                    r.fm_gp(j,i) = sqrt(flx(1,j)^2 + flx(2,j)^2);
                    
                    % Components of principal fluxes in X-Y directions
                    r.fmx_gp(npts) = r.fxx_gp(j,i);
                    r.fmy_gp(npts) = r.fyy_gp(j,i);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Extrapolate Gauss point results of derived variables to element
        % node results.
        % The nodal results are computed by extrapolation of Gauss point
        % results using the TGN matrix.
        % In thermal analysis, derived variables are heat fluxes.
        function elemDerivedVarExtrap(~,mdl)
            r = mdl.res;
            
            for i = 1:mdl.nel
                TGN = mdl.elems(i).TGN;
                nen = mdl.elems(i).shape.nen;
                ngp = r.ngp(i);
                
                r.fxx_elemextrap(1:nen,i) = TGN * r.fxx_gp(1:ngp,i);
                r.fyy_elemextrap(1:nen,i) = TGN * r.fyy_gp(1:ngp,i);
                r.fm_elemextrap(1:nen,i)  = TGN * r.fm_gp(1:ngp,i);
            end
        end
        
        %------------------------------------------------------------------
        % Smooth element node results of derived variables to global
        % node results.
        % The nodal global nodal results are computed by averaging values
        % of element extrapolated nodal results of all elements adjacent
        % to each node.
        % In thermal analysis, derived variables are heat fluxes.
        function nodeDerivedVarExtrap(~,mdl)
            r = mdl.res;
            
            % Sum contributions of extrapolated node results from connected elements
            for i = 1:mdl.nel
                for j = 1:mdl.elems(i).shape.nen
                    n = mdl.elems(i).shape.nodes(j).id;
                    r.fxx_nodeextrap(n) = r.fxx_nodeextrap(n) + r.fxx_elemextrap(j,i);
                    r.fyy_nodeextrap(n) = r.fyy_nodeextrap(n) + r.fyy_elemextrap(j,i);
                    r.fm_nodeextrap(n)  = r.fm_nodeextrap(n)  + r.fm_elemextrap(j,i);
                end
            end
            
            % Average nodal values by the number of connected elements
            for i = 1:mdl.nnp
                nadjelems = length(mdl.nodes(i).elems);
                r.fxx_nodeextrap(i) = r.fxx_nodeextrap(i) / nadjelems;
                r.fyy_nodeextrap(i) = r.fyy_nodeextrap(i) / nadjelems;
                r.fm_nodeextrap(i)  = r.fm_nodeextrap(i)  / nadjelems;
            end
        end
    end
end