%% Anm_Axisymmetric Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anm.html Anm: analysis model super-class> to deal
% with axisymmetric models.
%
classdef Anm_Axisymmetric < fem.Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm_Axisymmetric()
            this = this@fem.Anm(fem.Anm.AXISYMMETRIC,2);
            
            this.DISPL_X  = true;    % Displacement X
            this.DISPL_Y  = true;    % Displacement Y
            this.SIGMA_XX = true;    % Normal stress XX
            this.SIGMA_YY = true;    % Normal stress YY
            this.SIGMA_ZZ = true;    % Normal stress ZZ
            this.TAU_XY   = true;    % Shear stress XY
            this.SIGMA_1  = true;    % Principal stress 1
            this.SIGMA_2  = true;    % Principal stress 2
            this.TAU_MAX  = true;    % Maximum shear stress
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
                if (mdl.nodes(i).ebc(1) == 1)
                    mdl.neqc = mdl.neqc + 1;
                    mdl.ID(1,i) = 1;
                end

                % Check for fixed translation in global Y direction
                if (mdl.nodes(i).ebc(2) == 1)
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
                if (~isempty(mdl.nodes(i).prescDispl))
                    % Add prescribed displacement in global X direction
                    id = mdl.ID(1,i);
                    if (id > mdl.neqf)
                        D(id) = mdl.nodes(i).prescDispl(1);
                    end
                    
                    % Add prescribed displacement in global Y direction
                    id = mdl.ID(2,i);
                    if (id > mdl.neqf)
                        D(id) = mdl.nodes(i).prescDispl(2);
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
            
            C = e * [ 1-v  v    v    0;
                      v    1-v  v    0;
                      v    v    1-v  0;
                      0    0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble strain-displacement matrix at a given position of
        % an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        %  r,s: parametric coordinates
        function B = Bmtx(this,elem,GradNcar,r,s)                
            % Geometry map functions matrix evaluated at this point
            M = elem.shape.Mmtx(r,s);

            % Displacement shape functions matrix evaluated at this point
            N = elem.shape.Nmtx(r,s);
            
            % Location of evaluation point (X coordinate is the radius in axisymm.)
            p = M * elem.carCoord;
            radius = p(1);
            
            % Assemble strain-displacement matrix
            B = zeros(4,elem.shape.nen*this.ndof);
            
            for i = 1:elem.shape.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0.0;
                B(2,2*i-1) = N(i)/radius;     B(2,2*i) = 0.0;
                B(3,2*i-1) = 0.0;             B(3,2*i) = GradNcar(2,i);
                B(4,2*i-1) = GradNcar(2,i);   B(4,2*i) = GradNcar(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Returns the ridigity coefficient according to analysis model
        % at a given position in parametric coordinates of an element.
        % For axi-symmetric analysis, the rigidity coefficient is the
        % radius at the given point.
        function coeff = rigidityCoeff(~,elem,r,s)
            % Geometry map functions matrix evaluated at this point
            M = elem.shape.Mmtx(r,s);
            
            % Location of evaluation point (X coordinate is the radius in axisymm.)
            p = M * elem.carCoord;
            coeff = p(1);
        end
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        % Input:
        %  C:   constituive matrix
        %  B:   strain-displacement matrix
        %  d:   generalized displacements for all d.o.f.'s of element
        function str = pointStress(~,C,B,d)
            % Compute point stress components
            str_raw = C * B * d;
            
            % Skip tangential stress component
            str(1) = str_raw(1);
            str(2) = str_raw(3);
            str(3) = str_raw(4);
        end
    end
end