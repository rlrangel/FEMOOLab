%% Anm Class (Analysis Model)
%
%% Description
%
% This is an abstract super-class that generically specifies an analysis 
% model in the FEMOOLab program.
%
% Essentially, this super-class declares abstract methods that define
% the particular behavior of an analysis model. These abstract methods
% are the functions that should be implemented in a derived sub-class
% that deals with specific types of analysis models.
%
%% Subclasses
%
% * <anm_planestress.html Anm_PlaneStress: plane stress analysis model subclass>
% * <anm_planestrain.html Anm_PlaneStrain: plane strain analysis model subclass>
% * <anm_planeconduction.html Anm_PlaneConduction: plane conduction analysis model subclass>
% * <anm_axisymmetric.html Anm_Axisymmetric: axisymmetric analysis model subclass>
%
%% Class definition
%
classdef Anm < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of analysis model
        GENERIC          = int32(0);
        PLANE_STRESS     = int32(1);
        PLANE_STRAIN     = int32(2);
        PLANE_CONDUCTION = int32(3);
        AXISYMMETRIC     = int32(4);
    end
    
    %% Flags for types of responses
    properties (SetAccess = protected, GetAccess = public)
        % Elasticity properties
        DISPL_X     = logical(false);    % Displacement X
        DISPL_Y     = logical(false);    % Displacement Y
        DISPL_Z     = logical(false);    % Displacement Z
        ROTAT_X     = logical(false);    % Rotation about X
        ROTAT_Y     = logical(false);    % Rotation about Y
        ROTAT_Z     = logical(false);    % Rotation about Z
        SIGMA_XX    = logical(false);    % Normal stress XX
        SIGMA_YY    = logical(false);    % Normal stress YY
        SIGMA_ZZ    = logical(false);    % Normal stress ZZ
        TAU_XY      = logical(false);    % Shear stress XY
        TAU_XZ      = logical(false);    % Shear stress XZ
        TAU_YZ      = logical(false);    % Shear stress YZ
        SIGMA_1     = logical(false);    % Principal stress 1
        SIGMA_2     = logical(false);    % Principal stress 2
        SIGMA_3     = logical(false);    % Principal stress 3
        TAU_MAX     = logical(false);    % Maximum shear stress
        MOMENT_XX   = logical(false);    % Moment XX
        MOMENT_YY   = logical(false);    % Moment YY
        MOMENT_XY   = logical(false);    % Torsion moment XY
        SHEAR_XZ    = logical(false);    % Shear force XZ
        SHEAR_YZ    = logical(false);    % Shear force YZ
        MOMENT_1    = logical(false);    % Principal moment 1
        MOMENT_2    = logical(false);    % Principal moment 2
        TORSION_MAX = logical(false);    % Maximum torsion
        
        % Thermal properties
        TEMPERATURE = logical(false);    % Temperature
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;  % flag for type of analysis model
        ndof int32 = int32.empty;  % number of degrees-of-freedom per node
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm(type,ndof)
            this.type = type;
            this.ndof = ndof;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Initialize global d.o.f. numbering matrix with ones and zeros,
        % and count total number of equations of free and fixed  d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(j,i) = 0, d.o.f. j of node i is free.
        %  if ID(j,i) = 1, d.o.f. j of node i is fixed.
        setupDOFNum(this,mdl);
        
        %------------------------------------------------------------------
        % Assemble material constitutive matrix of a given element.
        C = Cmtx(this,elem);
        
        %------------------------------------------------------------------
        % Assemble strain matrix at a given position in parametric
        % coordinates of an element.
        B = Bmtx(this,elem,GradNcar,r,s);
        
        %------------------------------------------------------------------
        % Returns the ridigity coefficient at a given position in
        % parametric coordinates of an element.
        coeff = rigidityCoeff(this,elem,r,s);
        
        %------------------------------------------------------------------
        % Assemble global stiffness matrix.
        K = gblStiffMtx(this,mdl);
        
        %------------------------------------------------------------------
        % Add point force contributions to global forcing vector.
        F = addPointForce(this,mdl,F);
        
        %------------------------------------------------------------------
        % Add equivalent nodal force contributions to global forcing vector.
        F = addEquivForce(this,mdl,F);
        
        %------------------------------------------------------------------
        % Add essencial boundary conditions (prescribed values of state
        % variables) to global vector of state variables.
        U = addEBC(this,mdl,U);
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        str = pointStress(this,C,B,d);
    end
end