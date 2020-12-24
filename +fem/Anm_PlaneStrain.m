%% Anm_PlaneStrain Class (Plane Strain Model)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with plane strain models in an elasticity analysis.
%
%% Class definition
%
classdef Anm_PlaneStrain < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_PlaneStrain()
            this = this@fem.Anm(fem.Anm.PLANE_STRAIN,2);
            
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
        % Initialize global d.o.f. numbering matrix with ones and zeros,
        % and count total number of equations of free and fixed  d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(j,i) = 0, d.o.f. j of node i is free.
        %  if ID(j,i) = 1, d.o.f. j of node i is fixed.
        % Input arguments:
        %  mdl: handle to an object of the Model class
        function setupDOFNum(this,mdl)
            % Dimension global d.o.f. numbering matrix
            mdl.ID = zeros(this.ndof,mdl.nnp);
            
            % Initialize number of fixed d.o.f.'s
            mdl.neqc = 0;
            
            % Count number of fixed d.o.f.'s and setup ID matrix
            for i = 1:mdl.nnp
                 % Check for fixed translation in global X direction
                if (mdl.nodes(i).fixDispl(1))
                    mdl.neqc = mdl.neqc + 1;
                    mdl.ID(1,i) = 1;
                end

                % Check for fixed translation in global Y direction
                if (mdl.nodes(i).fixDispl(2))
                    mdl.neqc = mdl.neqc + 1;
                    mdl.ID(2,i) = 1;
                end
            end
            
            % Compute total number of free d.o.f.
            mdl.neqf = mdl.neq - mdl.neqc;
        end
        
        %------------------------------------------------------------------
        % Add point loads to global forcing vector, including the
        % components that correspond to fixed d.o.f.'s.
        % Input arguments:
        %  mdl: handle to an object of the Model class
        %  F: global forcing vector
        % Output arguments:
        %  F: global forcing vector
        function F = addPointLoad(~,mdl,F)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).load))
                     % Add applied force in global X direction
                    id = mdl.ID(1,i);
                    F(id) = mdl.nodes(i).load(1);
                    
                    % Add applied force in global Y direction
                    id = mdl.ID(2,i);
                    F(id) = mdl.nodes(i).load(2);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Adds prescribed displacements (known support settlement values)
        % to global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of global displacement vector that corresponds to a free d.o.f.
        % Input arguments:
        %  mdl: handle to an object of the Model class
        %  D: global displacement vector
        % Output arguments:
        %  D: global displacement vector
        function D = addPrescDispl(~,mdl,D)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).ebcDispl))
                    % Add prescribed displacement in global X direction
                    id = mdl.ID(1,i);
                    if (id > mdl.neqf)
                        D(id) = mdl.nodes(i).ebcDispl(1);
                    end
                    
                    % Add prescribed displacement in global Y direction
                    id = mdl.ID(2,i);
                    if (id > mdl.neqf)
                        D(id) = mdl.nodes(i).ebcDispl(2);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assemble material constitutive matrix for a given element.
        function C = Cmtx(~,elem)
            E = elem.mat.E;
            v = elem.mat.v;
            e = E/((1+v)*(1-(2*v)));
            
            C = e * [ 1-v  v    0;
                      v    1-v  0;
                      0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble strain-displacement matrix at a given position of
        % an element.
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
        % Returns the ridigity coefficient according to analysis model
        % at a given position in parametric coordinates of an element.
        % For plane strain analysis, the rigidity coefficient is one (1).
        function coeff = rigidityCoeff(~,~,~,~)
            coeff = 1;
        end
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        % Input:
        %  C:   constituive matrix
        %  B:   strain-displacement matrix
        %  d:   generalized displacements for all d.o.f.'s of element
        function str = pointStress(~,C,B,d)
            % In plane strain, raw stress vector is the target one
            str = C * B * d;
        end
    end
end