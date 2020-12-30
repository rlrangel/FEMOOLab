%% Anm_PlaneConduction Class (Plane Heat Conduction Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with plane heat conduction models in a thermal analysis.
%
%% Class definition
%
classdef Anm_PlaneConduction < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_PlaneConduction()
            this = this@fem.Anm(fem.Anm.PLANE_CONDUCTION,1);
            
            % Types of response
            this.TEMPERATURE = true;  % Temperature
            this.FLUX_XX     = true;  % flux XX
            this.FLUX_YY     = true;  % flux YY
            this.FLUX_PRC    = true;  % Principal flux
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anm
    methods
        %------------------------------------------------------------------
        % Initialize global d.o.f. numbering matrix with ones and zeros,
        % and count total number of equations of free and fixed  d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(j,i) = 0, d.o.f. j of node i is free.
        %  if ID(j,i) = 1, d.o.f. j of node i is fixed.
        function setupDOFNum(this,mdl)
            % Dimension global d.o.f. numbering matrix
            mdl.ID = zeros(this.ndof,mdl.nnp);
            
            % Initialize number of fixed d.o.f.'s
            mdl.neqc = 0;
            
            % Count number of fixed d.o.f.'s and setup ID matrix
            for i = 1:mdl.nnp
                % Check for fixed temperature
                if (mdl.nodes(i).fixTemp)
                    mdl.neqc = mdl.neqc + 1;
                    mdl.ID(1,i) = 1;
                end
            end
            
            % Compute total number of free d.o.f.
            mdl.neqf = mdl.neq - mdl.neqc;
        end
        
        %------------------------------------------------------------------
        % Assemble material constitutive matrix for a given element.
        function C = Cmtx(~,elem)
            k = elem.mat.k;
            
            C = [ k   0;
                  0   k ];
        end
        
        %------------------------------------------------------------------
        % Assemble strain matrix at a given position in parametric
        % coordinates of an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(this,elem,GradNcar,~,~)
            B = zeros(2,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,i) = GradNcar(1,i);
                B(2,i) = GradNcar(2,i);
            end
        end
        
        %------------------------------------------------------------------
        % Returns the ridigity coefficient at a given position in
        % parametric coordinates of an element.
        function coeff = rigidityCoeff(~,elem,~,~)
            coeff = elem.thk;
        end
        
        %------------------------------------------------------------------
        % Returns the mass coefficient.
        function coeff = massCoeff(~,elem)
            coeff = elem.rho * elem.cp;
        end
        
        %------------------------------------------------------------------
        % Assemble global stiffness matrix.
        function K = gblStiffMtx(~,mdl)
            % Initialize global stiffness matrix
            K = zeros(mdl.neq,mdl.neq);
            
            % Get element stiffness matrices and assemble global matrix
            for i = 1:mdl.nel
                gle = mdl.elems(i).gle;
                
                % Conventional stiffness
                ke = mdl.elems(i).stiffMtx();
                
                % Convection stiffness
                if (~isempty(mdl.elems(i).lineConvec))
                    ke = ke + mdl.elems(i).stiffConvecMtx();
                end
                
                K(gle,gle) = K(gle,gle) + ke;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to first time derivative of state
        % variables (capacity matrix in thermal analysis).
        function C = gblVelMtx(~,mdl)
            % Initialize global capacity matrix
            C = zeros(mdl.neq,mdl.neq);
            
            % Get element capacity matrices and assemble global matrix
            for i = 1:mdl.nel
                gle = mdl.elems(i).gle;
                ce = mdl.elems(i).massMtx();
                C(gle,gle) = C(gle,gle) + ce;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to second time derivative of state
        % variables ("acceleration" matrix).
        function M = gblAccelMtx(~,~)
            % NOT IMPLEMENTED
            M = zeros(mdl.neq,mdl.neq);
        end
        
        %------------------------------------------------------------------
        % Assemble global initial conditions matrix.
        function IC = gblInitCondMtx(~,mdl)
            % Initialize initial conditions matrix
            IC = zeros(mdl.neqf,1);
            
            for i = 1:mdl.nnp
                id  = mdl.ID(1,i);
                
                % Apply initial temperature only to free d.o.f.'s
                if (id <= mdl.neqf)
                    % If initial temperature is not given, assume zero
                    if (~isempty(mdl.nodes(i).iniTemp))
                        IC(id,1) = mdl.nodes(i).iniTemp;
                    else
                        IC(id,1) = 0;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add point force contributions to global forcing vector.
        function F = addPointForce(~,mdl,F)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).flux))
                    % Add load to global load vector
                    id  = mdl.ID(1,i);
                    F(id) = F(id) + mdl.nodes(i).flux;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add equivalent nodal force contributions to global forcing vector.
        function F = addEquivForce(~,mdl,F)
            for i = 1:mdl.nel
                gle = mdl.elems(i).gle;
                
                % Get element equivalent nodal flux vectors and assemble global vector
                if (~isempty(mdl.elems(i).lineFlux))
                    fline = mdl.elems(i).edgeEquivForceVct(mdl.elems(i).lineFlux);
                    F(gle) = F(gle) + fline;
                end
                
                if (~isempty(mdl.elems(i).lineConvec))
                    % New matrix whose last column is the product of
                    % convection coeff. and ambient temperature
                    % (columns 3 and 4 of property lineConvec)
                    lineConvec = mdl.elems(i).lineConvec(:,1:3);
                    lineConvec(:,3) = lineConvec(:,3) .* mdl.elems(i).lineConvec(:,4);
                    
                    fconv = mdl.elems(i).edgeEquivForceVct(lineConvec);
                    F(gle) = F(gle) + fconv;
                end
                
                if (~isempty(mdl.elems(i).domainFlux))
                    fdom = mdl.elems(i).domainEquivForceVct(mdl.elems(i).domainFlux);
                    F(gle) = F(gle) + fdom;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add essencial boundary conditions (prescribed values of state
        % variables) to global vector of state variables.
        % Avoid storing a prescribed value in a position of global vector
        % of state variables that corresponds to a free d.o.f.
        function U = addEBC(~,mdl,U)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).ebcTemp))
                    % Add prescribed temperature to global vector
                    id = mdl.ID(1,i);
                    if (id > mdl.neqf)
                        U(id) = mdl.nodes(i).ebcTemp;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        % Input:
        %  C:   constituive matrix
        %  B:   strain-displacement matrix
        %  d:   generalized displacements for all d.o.f.'s of element
        function str = pointStress(~,C,B,d)
            % In plane stress, raw stress vector is the target one
            str = -C * B * d;
        end
    end
end