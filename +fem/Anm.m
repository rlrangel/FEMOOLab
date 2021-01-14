%% Anm Class (Analysis Model)
%
%% Description
%
% This is an abstract super-class that generically specifies an analysis 
% model in the FEMOOLab program.
% Essentially, this super-class declares abstract methods that define
% the general behavior of analysis models. These abstract methods
% are the functions that should be implemented in a derived sub-class
% that deals with specific types of analysis models.
%
%% Subclasses
%
% * <anm_planestress.html Anm_PlaneStress: plane stress analysis model subclass>
% * <anm_planestrain.html Anm_PlaneStrain: plane strain analysis model subclass>
% * <anm_axisymstress.html Anm_AxisymStress: axisymmetric stress analysis model subclass>
% * <anm_thickplate.html Anm_ThickPlate: thick plate analysis model subclass>
% * <anm_planeconduction.html Anm_PlaneConduction: plane heat conduction analysis model subclass>
% * <anm_axisymconduction.html Anm_AxisymConduction: axisymmetric heat conduction analysis model subclass>
% * <anm_planeconvdiff.html Anm_PlaneConvDiff: plane heat convection diffusion analysis model subclass>
%
%% Class definition
%
classdef Anm < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of physics
        STRUCTURAL = int32(1);
        THERMAL    = int32(2);
        
        % Types of analysis model
        PLANE_STRESS         = int32(1);
        PLANE_STRAIN         = int32(2);
        AXISYM_STRESS        = int32(3);
        THICK_PLATE          = int32(4);
        PLANE_CONDUCTION     = int32(5);
        AXISYM_CONDUCTION    = int32(6);
        CONVECTION_DIFFUSION = int32(7);
    end
    
    %% Flags for types of responses
    properties (SetAccess = protected, GetAccess = public)
        % Structural properties
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
        SHEAR_XZ    = logical(false);    % Shear force XZ
        SHEAR_YZ    = logical(false);    % Shear force YZ
        MOMENT_XX   = logical(false);    % Moment XX
        MOMENT_YY   = logical(false);    % Moment YY
        MOMENT_XY   = logical(false);    % Torsion moment XY
        MOMENT_1    = logical(false);    % Principal moment 1
        MOMENT_2    = logical(false);    % Principal moment 2
        TORSION_MAX = logical(false);    % Maximum torsion moment
        
        % Thermal properties
        TEMPERATURE = logical(false);    % Temperature
        FLUX_XX     = logical(false);    % Heat flux XX
        FLUX_YY     = logical(false);    % Heat flux YY
        FLUX_ZZ     = logical(false);    % Heat flux ZZ
        FLUX_MOD    = logical(false);    % Heat principal flux
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        phys int32 = int32.empty;  % flag for type of physics
        type int32 = int32.empty;  % flag for type of analysis model
        ndof int32 = int32.empty;  % number of d.o.f.'s per node
        ndvc int32 = int32.empty;  % number of derived variable components
        gla  int32 = int32.empty;  % gather vector (stores local displ. d.o.f. numbers of a node)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anm(phys,type,ndof,ndvc,gla)
            this.phys = phys;
            this.type = type;
            this.ndof = ndof;
            this.ndvc = ndvc;
            this.gla  = gla;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Assemble material constitutive matrix of a given element.
        C = Cmtx(this,elem);
        
        %------------------------------------------------------------------
        % Assemble strain matrix at a given position in parametric
        % coordinates of an element.
        B = Bmtx(this,elem,GradNcar,r,s);
        
        %------------------------------------------------------------------
        % Return the ridigity coefficient at a given position in
        % parametric coordinates of an element.
        coeff = rigidityCoeff(this,elem,r,s);
        
        %------------------------------------------------------------------
        % Return the mass coefficient.
        coeff = massCoeff(this,elem);
        
        %------------------------------------------------------------------
        % Assemble global stiffness matrix.
        K = gblStiffMtx(this,mdl);
        
        %------------------------------------------------------------------
        % Assemble global matrix related to 1st time derivative of d.o.f.'s.
        C = gblRate1Mtx(this,mdl);
        
        %------------------------------------------------------------------
        % Assemble global matrix related to 2nd time derivative of d.o.f.'s.
        M = gblRate2Mtx(this,mdl);
        
        %------------------------------------------------------------------
        % Modify system arrays to include stabilization components for the
        % convective term.
        [K,F] = stabConvec(this,mdl,K,F);
        
        %------------------------------------------------------------------
        % Compute derived variable components at a given point of an element.
        dvar = pointDerivedVar(this,C,B,d);
        
        %------------------------------------------------------------------
        % Compute derived variables and the principal values and directions
        % at Gauss points of all elements.
        gaussDerivedVar(this,mdl);
        
        %------------------------------------------------------------------
        % Extrapolate Gauss point results of derived variables to element
        % node results.
        % The nodal results are computed by extrapolation of Gauss point
        % results using the TGN matrix.
        elemDerivedVarExtrap(this,mdl);
        
        %------------------------------------------------------------------
        % Smooth element node results of derived variables to global
        % node results.
        % The nodal global nodal results are computed by averaging values
        % of element extrapolated nodal results of all elements adjacent
        % to each node.
        nodeDerivedVarExtrap(this,mdl);
    end
end