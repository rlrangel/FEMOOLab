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
        % Add point force contributions to global forcing vector,
        % including the components that correspond to fixed d.o.f.'s.
        function F = addPointForce(~,~,F)
            % Point flux is not considered
            return;
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
            str = C * B * d;
        end
    end
end