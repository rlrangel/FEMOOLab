%% Anm_PlaneStress Class (Plane Stress Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with plane stress models in a structural analysis.
%
%% Class definition
%
classdef Anm_PlaneStress < fem.Anm
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        gla int32 = int32.empty;  % gather vector (stores local displ. d.o.f. numbers of a node)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_PlaneStress()
            this = this@fem.Anm(fem.Anm.PLANE_STRESS,2);
            this.gla = [1 2];
            
            % Types of response
            this.DISPL_X  = true;  % Displacement X
            this.DISPL_Y  = true;  % Displacement Y
            this.SIGMA_XX = true;  % Normal stress XX
            this.SIGMA_YY = true;  % Normal stress YY
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
                for j = 1:this.ndof
                    % Get d.o.f number
                    dof = this.gla(j);
                    
                    % Check for fixed d.o.f
                    if (mdl.nodes(i).fixDispl(dof))
                        mdl.neqc = mdl.neqc + 1;
                        mdl.ID(j,i) = 1;
                    end
                end
            end
            
            % Compute total number of free d.o.f.
            mdl.neqf = mdl.neq - mdl.neqc;
        end
        
        %------------------------------------------------------------------
        % Assemble material constitutive matrix for a given element.
        function C = Cmtx(~,elem)
            E = elem.mat.E;
            v = elem.mat.v;
            e = E/(1-(v^2));
            
            C = e * [ 1    v    0;
                      v    1    0;
                      0    0    (1-v)/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble strain matrix at a given position in parametric
        % coordinates of an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(this,elem,GradNcar,~,~)
            B = zeros(3,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0;
                B(2,2*i-1) = 0;               B(2,2*i) = GradNcar(2,i);
                B(3,2*i-1) = GradNcar(2,i);   B(3,2*i) = GradNcar(1,i);
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
            coeff = elem.rho;
        end
        
        %------------------------------------------------------------------
        % Assemble global stiffness matrix.
        function K = gblStiffMtx(~,mdl)
            % Initialize global stiffness matrix
            K = zeros(mdl.neq,mdl.neq);
            
            % Get element stiffness matrices and assemble global matrix
            for i = 1:mdl.nel
                gle = mdl.elems(i).gle;
                ke = mdl.elems(i).stiffMtx();
                K(gle,gle) = K(gle,gle) + ke;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to first time derivative of state
        % variables (damping matrix in structural analysis).
        function C = gblVelMtx(~,~)
            % NOT IMPLEMENTED
            C = zeros(mdl.neq,mdl.neq);
        end
        
        %------------------------------------------------------------------
        % Assemble global matrix related to second time derivative of state
        % variables (mass matrix in structural analysis).
        function M = gblAccelMtx(~,mdl)
            % Initialize global mass matrix
            M = zeros(mdl.neq,mdl.neq);
            
            % Get element mass matrices and assemble global matrix
            for i = 1:mdl.nel
                gle = mdl.elems(i).gle;
                me = mdl.elems(i).massMtx();
                M(gle,gle) = M(gle,gle) + me;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global initial conditions matrix.
        function IC = gblInitCondMtx(~,mdl)
            % Initialize initial conditions matrix
            IC = zeros(mdl.neqf,2);
            
            for i = 1:mdl.nnp
                for j = 1:mdl.anm.ndof
                    % Get d.o.f numbers
                    id  = mdl.ID(j,i);
                    dof = mdl.anm.gla(j);
                    
                    % Apply initial conditions only to free d.o.f.'s
                    if (id <= mdl.neqf)
                        if (~isempty(mdl.nodes(i).iniDispl) ||...
                            ~isempty(mdl.nodes(i).iniVeloc))
                            IC(id,1) = mdl.nodes(i).iniDispl(dof);
                            IC(id,2) = mdl.nodes(i).iniVeloc(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add point force contributions to global forcing vector.
        function F = addPointForce(this,mdl,F)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).load))
                    for j = 1:this.ndof
                        % Get d.o.f numbers
                        id  = mdl.ID(j,i);
                        dof = this.gla(j);
                        
                        % Add load to reference load vector
                        F(id) = F(id) + mdl.nodes(i).load(dof);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add equivalent nodal force contributions to global forcing vector.
        function F = addEquivForce(~,mdl,F)
            for i = 1:mdl.nel
                gle = mdl.elems(i).gle;
                
                % Get element equivalent nodal load vectors and assemble global vector
                if (~isempty(mdl.elems(i).lineLoad))
                    fline = mdl.elems(i).edgeEquivForceVct(mdl.elems(i).lineLoad);
                    F(gle) = F(gle) + fline;
                end
                
                if (~isempty(mdl.elems(i).domainLoad))
                    fdom = mdl.elems(i).domainEquivForceVct(mdl.elems(i).domainLoad);
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
                if (~isempty(mdl.nodes(i).ebcDispl))
                    for j = 1:this.ndof
                        % Get d.o.f numbers
                        id  = mdl.ID(j,i);
                        dof = this.gla(j);
                        
                        % Add prescribed displacement to global vector
                        if (id > mdl.neqf)
                            U(id) = mdl.nodes(i).ebcDispl(dof);
                        end
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
            str = C * B * d;
        end
    end
end