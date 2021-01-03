%% Result Class
%
%% Description
%
% This class defines a result object in the FEMOOLab program.
% A Result object is responsible for storing the analysis results and
% perform some operations on them.
%
classdef Result < handle
    %% Flags for response plotting options
    properties (SetAccess = public, GetAccess = public)
        % General
        eid    = logical(false);  % element numbers
        nid    = logical(false);  % node numbers
        gid    = logical(false);  % gauss numbers
        smooth = logical(false);  % Smooth element results at common nodes
        tol    = double(1e-5);    % tolerance for cleaning small result values and differences
        
        % Structural analysis contours
        scl    = double(0.0);     % Scale factor for deformed mesh
        deform = logical(false);  % deformed mesh
        dx     = logical(false);  % displacements in X direction
        dy     = logical(false);  % displacements in Y direction
        dz     = logical(false);  % displacements in Z directiont
        rx     = logical(false);  % rotations about X axis
        ry     = logical(false);  % rotations about Y axis
        rz     = logical(false);  % rotations about Z axis
        sxx    = logical(false);  % normal stresses in X direction
        syy    = logical(false);  % normal stresses in Y direction
        szz    = logical(false);  % normal stresses in Z direction
        txy    = logical(false);  % XY shear stresses
        txz    = logical(false);  % XZ shear stresses
        tyz    = logical(false);  % YZ shear stresses
        s1     = logical(false);  % principal stresses 1
        s2     = logical(false);  % principal stresses 2
        s3     = logical(false);  % principal stresses 3
        taumax = logical(false);  % maximum shear stresses
        mxx    = logical(false);  % moment about X direction
        myy    = logical(false);  % moment about Y direction
        mxy    = logical(false);  % torsion moment
        qxz    = logical(false);  % XZ shear force
        qyz    = logical(false);  % YZ shear force
        m1     = logical(false);  % principal moment 1
        m2     = logical(false);  % principal moment 2
        tormax = logical(false);  % maximum torsion
        
        % Thermal analysis contours
        temp   = logical(false);  % temperature field
        fxx    = logical(false);  % heat fluxes in X direction
        fyy    = logical(false);  % heat fluxes in Y direction
        fzz    = logical(false);  % heat fluxes in Z direction
        fm     = logical(false);  % heat fluxes module
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        maxGPts int32  = int32.empty;  % max. number of Gauss points in all elements
        maxNen  int32  = int32.empty;  % max. number of nodal points in all elements
        
        % Results from equilibrium system
        U     double = double.empty;   % global vector of state variables
        Ut    double = double.empty;   % global vector of state variable 1st time derivatives
        Utt   double = double.empty;   % global vector of state variable 2nd time derivatives
        steps int32  = int32.empty;    % number of performed steps
        times double = double.empty;   % vector of time values of each step
        
        % Gauss point results
        ngp                 int32  = int32.empty;  % vector of number of element gauss pts
        x_gp                double = double.empty; % vector of gauss point x coordinates
        y_gp                double = double.empty; % vector of gauss point y coordinates
        z_gp                double = double.empty; % vector of gauss point Z coordinates
        
        sxx_gp              double = double.empty; % sigma x gauss points stress array
        sxx_gp_min          double = double.empty; % sigma x gauss points min. stress
        sxx_gp_max          double = double.empty; % sigma x gauss points max. stress
        syy_gp              double = double.empty; % sigma y gauss points stress array
        syy_gp_min          double = double.empty; % sigma y gauss points min. stress
        syy_gp_max          double = double.empty; % sigma y gauss points max. stress
        szz_gp              double = double.empty; % sigma Z gauss points stress array
        szz_gp_min          double = double.empty; % sigma Z gauss points min. stress
        szz_gp_max          double = double.empty; % sigma Z gauss points max. stress
        txy_gp              double = double.empty; % tau xy gauss points stress array
        txy_gp_min          double = double.empty; % tau xy gauss points min. stress
        txy_gp_max          double = double.empty; % tau xy gauss points max. stress
        txz_gp              double = double.empty; % tau xz gauss points stress array
        txz_gp_min          double = double.empty; % tau xz gauss points min. stress
        txz_gp_max          double = double.empty; % tau xz gauss points max. stress
        tyz_gp              double = double.empty; % tau yz gauss points stress array
        tyz_gp_min          double = double.empty; % tau yz gauss points min. stress
        tyz_gp_max          double = double.empty; % tau yz gauss points max. stress
        s1_gp               double = double.empty; % sigma 1 gauss points stress array
        s1_gp_min           double = double.empty; % sigma 1 gauss points min. stress
        s1_gp_max           double = double.empty; % sigma 1 gauss points max. stress
        s2_gp               double = double.empty; % sigma 2 gauss points stress array
        s2_gp_min           double = double.empty; % sigma 2 gauss points min. stress
        s2_gp_max           double = double.empty; % sigma 2 gauss points max. stress
        s3_gp               double = double.empty; % sigma 3 gauss points stress array
        s3_gp_min           double = double.empty; % sigma 3 gauss points min. stress
        s3_gp_max           double = double.empty; % sigma 3 gauss points max. stress
        tmax_gp             double = double.empty; % tau max. gauss points stress array
        tmax_gp_min         double = double.empty; % tau max. gauss points min. stress
        tmax_gp_max         double = double.empty; % tau max. gauss points max. stress
        s1x_gp              double = double.empty; % vector of s1 vector gauss x components
        s1y_gp              double = double.empty; % vector of s1 vector gauss y components
        s1z_gp              double = double.empty; % vector of s1 vector gauss z components
        s2x_gp              double = double.empty; % vector of s2 vector gauss x components
        s2y_gp              double = double.empty; % vector of s2 vector gauss y components
        s2z_gp              double = double.empty; % vector of s2 vector gauss y components
        s3x_gp              double = double.empty; % vector of s3 vector gauss x components
        s3y_gp              double = double.empty; % vector of s3 vector gauss y components
        s3z_gp              double = double.empty; % vector of s3 vector gauss z components
        
        fxx_gp              double = double.empty; % heat flux x gauss points array
        fxx_gp_min          double = double.empty; % heat flux x gauss points min. value
        fxx_gp_max          double = double.empty; % heat flux x gauss points max. value
        fyy_gp              double = double.empty; % heat flux y gauss points array
        fyy_gp_min          double = double.empty; % heat flux y gauss points min. value
        fyy_gp_max          double = double.empty; % heat flux y gauss points max. value
        fzz_gp              double = double.empty; % heat flux z gauss points array
        fzz_gp_min          double = double.empty; % heat flux z gauss points min. value
        fzz_gp_max          double = double.empty; % heat flux z gauss points max. value
        fm_gp               double = double.empty; % heat flux module gauss points array
        fm_gp_min           double = double.empty; % heat flux module gauss points min. value
        fm_gp_max           double = double.empty; % heat flux module gauss points max. value
        fmx_gp              double = double.empty; % vector of gauss heat flux module x components
        fmy_gp              double = double.empty; % vector of gauss heat flux module y components
        fmz_gp              double = double.empty; % vector of gauss heat flux module z components
        
        % Element node results
        sxx_elemextrap      double = double.empty; % sigma x element node extrap. stress array
        sxx_elemextrap_min  double = double.empty; % sigma x element node extrap. min. stress
        sxx_elemextrap_max  double = double.empty; % sigma x element node extrap. max. stress
        syy_elemextrap      double = double.empty; % sigma y element node extrap. stress array
        syy_elemextrap_min  double = double.empty; % sigma y element node extrap. min. stress
        syy_elemextrap_max  double = double.empty; % sigma y element node extrap. max. stress
        szz_elemextrap      double = double.empty; % sigma z element node extrap. stress array
        szz_elemextrap_min  double = double.empty; % sigma z element node extrap. min. stress
        szz_elemextrap_max  double = double.empty; % sigma z element node extrap. max. stress
        txy_elemextrap      double = double.empty; % tau xy element node extrap. stress array
        txy_elemextrap_min  double = double.empty; % tau xy element node extrap. min. stress
        txy_elemextrap_max  double = double.empty; % tau xy element node extrap. max. stress
        txz_elemextrap      double = double.empty; % tau xz element node extrap. stress array
        txz_elemextrap_min  double = double.empty; % tau xz element node extrap. min. stress
        txz_elemextrap_max  double = double.empty; % tau xz element node extrap. max. stress
        tyz_elemextrap      double = double.empty; % tau yz element node extrap. stress array
        tyz_elemextrap_min  double = double.empty; % tau yz element node extrap. min. stress
        tyz_elemextrap_max  double = double.empty; % tau yz element node extrap. max. stress
        s1_elemextrap       double = double.empty; % sigma 1 element node extrap. stress array
        s1_elemextrap_min   double = double.empty; % sigma 1 element node extrap. min. stress
        s1_elemextrap_max   double = double.empty; % sigma 1 element node extrap. max. stress
        s2_elemextrap       double = double.empty; % sigma 2 element node extrap. stress array
        s2_elemextrap_min   double = double.empty; % sigma 2 element node extrap. min. stress
        s2_elemextrap_max   double = double.empty; % sigma 2 element node extrap. max. stress
        s3_elemextrap       double = double.empty; % sigma 3 element node extrap. stress array
        s3_elemextrap_min   double = double.empty; % sigma 3 element node extrap. min. stress
        s3_elemextrap_max   double = double.empty; % sigma 3 element node extrap. max. stress
        tmax_elemextrap     double = double.empty; % tau max. element node extrap. stress array
        tmax_elemextrap_min double = double.empty; % tau max. element node extrap. min stress
        tmax_elemextrap_max double = double.empty; % tau max. element node extrap. max stress
        
        fxx_elemextrap      double = double.empty; % flux x element node extrap. array
        fxx_elemextrap_min  double = double.empty; % flux x element node extrap. min. value
        fxx_elemextrap_max  double = double.empty; % flux x element node extrap. max. value
        fyy_elemextrap      double = double.empty; % flux y element node extrap. array
        fyy_elemextrap_min  double = double.empty; % flux y element node extrap. min. value
        fyy_elemextrap_max  double = double.empty; % flux y element node extrap. max. value
        fzz_elemextrap      double = double.empty; % flux z element node extrap. array
        fzz_elemextrap_min  double = double.empty; % flux z element node extrap. min. value
        fzz_elemextrap_max  double = double.empty; % flux z element node extrap. max. value
        fm_elemextrap       double = double.empty; % flux module element node extrap. array
        fm_elemextrap_min   double = double.empty; % flux module element node extrap. min. value
        fm_elemextrap_max   double = double.empty; % flux module element node extrap. max. value
        
        % Global node results
        sxx_nodeextrap      double = double.empty; % sigma x extrap. node smoothed stress array
        sxx_nodeextrap_min  double = double.empty; % sigma x extrap. node smoothed min. stress
        sxx_nodeextrap_max  double = double.empty; % sigma x extrap. node smoothed max. stress
        syy_nodeextrap      double = double.empty; % sigma y extrap. node smoothed stress array
        syy_nodeextrap_min  double = double.empty; % sigma y extrap. node smoothed min. stress
        syy_nodeextrap_max  double = double.empty; % sigma y extrap. node smoothed max. stress
        szz_nodeextrap      double = double.empty; % sigma z extrap. node smoothed stress array
        szz_nodeextrap_min  double = double.empty; % sigma z extrap. node smoothed min. stress
        szz_nodeextrap_max  double = double.empty; % sigma z extrap. node smoothed max. stress
        txy_nodeextrap      double = double.empty; % tau xy extrap. node smoothed stress array
        txy_nodeextrap_min  double = double.empty; % tau xy extrap. node smoothed min. stress
        txy_nodeextrap_max  double = double.empty; % tau xy extrap. node smoothed max. stress
        txz_nodeextrap      double = double.empty; % tau xz extrap. node smoothed stress array
        txz_nodeextrap_min  double = double.empty; % tau xz extrap. node smoothed min. stress
        txz_nodeextrap_max  double = double.empty; % tau xz extrap. node smoothed max. stress
        tyz_nodeextrap      double = double.empty; % tau yz extrap. node smoothed stress array
        tyz_nodeextrap_min  double = double.empty; % tau yz extrap. node smoothed min. stress
        tyz_nodeextrap_max  double = double.empty; % tau yz extrap. node smoothed max. stress
        s1_nodeextrap       double = double.empty; % sigma 1 extrap. node smoothed stress array
        s1_nodeextrap_min   double = double.empty; % sigma 1 extrap. node smoothed min. stress
        s1_nodeextrap_max   double = double.empty; % sigma 1 extrap. node smoothed max. stress
        s2_nodeextrap       double = double.empty; % sigma 2 extrap. node smoothed stress array
        s2_nodeextrap_min   double = double.empty; % sigma 2 extrap. node smoothed min. stress
        s2_nodeextrap_max   double = double.empty; % sigma 2 extrap. node smoothed max. stress
        s3_nodeextrap       double = double.empty; % sigma 3 extrap. node smoothed stress array
        s3_nodeextrap_min   double = double.empty; % sigma 3 extrap. node smoothed min. stress
        s3_nodeextrap_max   double = double.empty; % sigma 3 extrap. node smoothed max. stress
        tmax_nodeextrap     double = double.empty; % tau max. extrap. node smoothed stress array
        tmax_nodeextrap_min double = double.empty; % tau max. extrap. node smoothed min. stress
        tmax_nodeextrap_max double = double.empty; % tau max. extrap. node smoothed max. stress
        
        fxx_nodeextrap      double = double.empty; % flux x extrap. node smoothed array
        fxx_nodeextrap_min  double = double.empty; % flux x extrap. node smoothed min. value
        fxx_nodeextrap_max  double = double.empty; % flux x extrap. node smoothed max. value
        fyy_nodeextrap      double = double.empty; % flux y extrap. node smoothed array
        fyy_nodeextrap_min  double = double.empty; % flux y extrap. node smoothed min. value
        fyy_nodeextrap_max  double = double.empty; % flux y extrap. node smoothed max. value
        fzz_nodeextrap      double = double.empty; % flux z extrap. node smoothed array
        fzz_nodeextrap_min  double = double.empty; % flux z extrap. node smoothed min. value
        fzz_nodeextrap_max  double = double.empty; % flux z extrap. node smoothed max. value
        fm_nodeextrap       double = double.empty; % flux module extrap. node smoothed array
        fm_nodeextrap_min   double = double.empty; % flux module extrap. node smoothed min. value
        fm_nodeextrap_max   double = double.empty; % flux module extrap. node smoothed max. value
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Result()
            return;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Initialize arrays of pos-process results.
        function initPosResults(this,mdl)
            % Get maximum number of Gauss points and nodes of all elements
            this.maxGPts = mdl.maxGaussDerivedVarNpts();
            this.maxNen  = mdl.maxNumElemNodes();
            
            % Initialize arrays of active results
            this.ngp  = zeros(mdl.nel,1);
            this.x_gp = zeros(this.maxGPts*mdl.nel,1);
            this.y_gp = zeros(this.maxGPts*mdl.nel,1);
            this.z_gp = zeros(this.maxGPts*mdl.nel,1);
            
            if (this.sxx)
                this.sxx_gp         = zeros(this.maxGPts,mdl.nel);
                this.sxx_elemextrap = zeros(this.maxNen,mdl.nel);
                this.sxx_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.syy)
                this.syy_gp         = zeros(this.maxGPts,mdl.nel);
                this.syy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.syy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.szz)
                this.szz_gp         = zeros(this.maxGPts,mdl.nel);
                this.szz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.szz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.txy)
                this.txy_gp         = zeros(this.maxGPts,mdl.nel);
                this.txy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.txy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.txz)
                this.txz_gp         = zeros(this.maxGPts,mdl.nel);
                this.txz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.txz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.tyz)
                this.tyz_gp         = zeros(this.maxGPts,mdl.nel);
                this.tyz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.tyz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.s1)
                this.s1_gp         = zeros(this.maxGPts,mdl.nel);
                this.s1x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s1y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s1_elemextrap = zeros(this.maxNen,mdl.nel);
                this.s1_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.s2)
                this.s2_gp         = zeros(this.maxGPts,mdl.nel);
                this.s2x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s2y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s2_elemextrap = zeros(this.maxNen,mdl.nel);
                this.s2_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.s3)
                this.s3_gp         = zeros(this.maxGPts,mdl.nel);
                this.s3x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s3y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s3_elemextrap = zeros(this.maxNen,mdl.nel);
                this.s3_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.taumax)
                this.tmax_gp         = zeros(this.maxGPts,mdl.nel);
                this.tmax_elemextrap = zeros(this.maxNen,mdl.nel);
                this.tmax_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.fxx)
                this.fxx_gp         = zeros(this.maxGPts,mdl.nel);
                this.fxx_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fxx_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.fyy)
                this.fyy_gp         = zeros(this.maxGPts,mdl.nel);
                this.fyy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fyy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.fzz)
                this.fzz_gp         = zeros(this.maxGPts,mdl.nel);
                this.fzz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fzz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (this.fm)
                this.fm_gp         = zeros(this.maxGPts,mdl.nel);
                this.fmx_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.fmy_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.fm_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fm_nodeextrap = zeros(mdl.nnp,1);
            end
        end
        
        %------------------------------------------------------------------
        % Compute minimum and maximum values of obtained results.
        function setMinMaxValues(this)
            if (this.sxx)
                this.sxx_gp_min         = min(min(this.sxx_gp));
                this.sxx_gp_max         = max(max(this.sxx_gp));
                this.sxx_elemextrap_min = min(min(this.sxx_elemextrap));
                this.sxx_elemextrap_max = max(max(this.sxx_elemextrap));
                this.sxx_nodeextrap_min = min(min(this.sxx_nodeextrap));
                this.sxx_nodeextrap_max = max(max(this.sxx_nodeextrap));
            end
            if (this.syy)
                this.syy_gp_min         = min(min(this.syy_gp));
                this.syy_gp_max         = max(max(this.syy_gp));
                this.syy_elemextrap_min = min(min(this.syy_elemextrap));
                this.syy_elemextrap_max = max(max(this.syy_elemextrap));
                this.syy_nodeextrap_min = min(min(this.syy_nodeextrap));
                this.syy_nodeextrap_max = max(max(this.syy_nodeextrap));
            end
            if (this.szz)
                this.szz_gp_min         = min(min(this.szz_gp));
                this.szz_gp_max         = max(max(this.szz_gp));
                this.szz_elemextrap_min = min(min(this.szz_elemextrap));
                this.szz_elemextrap_max = max(max(this.szz_elemextrap));
                this.szz_nodeextrap_min = min(min(this.szz_nodeextrap));
                this.szz_nodeextrap_max = max(max(this.szz_nodeextrap));
            end
            if (this.txy)
                this.txy_gp_min         = min(min(this.txy_gp));
                this.txy_gp_max         = max(max(this.txy_gp));
                this.txy_elemextrap_min = min(min(this.txy_elemextrap));
                this.txy_elemextrap_max = max(max(this.txy_elemextrap));
                this.txy_nodeextrap_min = min(min(this.txy_nodeextrap));
                this.txy_nodeextrap_max = max(max(this.txy_nodeextrap));
            end
            if (this.txz)
                this.txz_gp_min         = min(min(this.txz_gp));
                this.txz_gp_max         = max(max(this.txz_gp));
                this.txz_elemextrap_min = min(min(this.txz_elemextrap));
                this.txz_elemextrap_max = max(max(this.txz_elemextrap));
                this.txz_nodeextrap_min = min(min(this.txz_nodeextrap));
                this.txz_nodeextrap_max = max(max(this.txz_nodeextrap));
            end
            if (this.tyz)
                this.tyz_gp_min         = min(min(this.tyz_gp));
                this.tyz_gp_max         = max(max(this.tyz_gp));
                this.tyz_elemextrap_min = min(min(this.tyz_elemextrap));
                this.tyz_elemextrap_max = max(max(this.tyz_elemextrap));
                this.tyz_nodeextrap_min = min(min(this.tyz_nodeextrap));
                this.tyz_nodeextrap_max = max(max(this.tyz_nodeextrap));
            end
            if (this.s1)
                this.s1_gp_min         = min(min(this.s1_gp));
                this.s1_gp_max         = max(max(this.s1_gp));
                this.s1_elemextrap_min = min(min(this.s1_elemextrap));
                this.s1_elemextrap_max = max(max(this.s1_elemextrap));
                this.s1_nodeextrap_min = min(min(this.s1_nodeextrap));
                this.s1_nodeextrap_max = max(max(this.s1_nodeextrap));
            end
            if (this.s2)
                this.s2_gp_min         = min(min(this.s2_gp));
                this.s2_gp_max         = max(max(this.s2_gp));
                this.s2_elemextrap_min = min(min(this.s2_elemextrap));
                this.s2_elemextrap_max = max(max(this.s2_elemextrap));
                this.s2_nodeextrap_min = min(min(this.s2_nodeextrap));
                this.s2_nodeextrap_max = max(max(this.s2_nodeextrap));
            end
            if (this.s3)
                this.s3_gp_min         = min(min(this.s3_gp));
                this.s3_gp_max         = max(max(this.s3_gp));
                this.s3_elemextrap_min = min(min(this.s3_elemextrap));
                this.s3_elemextrap_max = max(max(this.s3_elemextrap));
                this.s3_nodeextrap_min = min(min(this.s3_nodeextrap));
                this.s3_nodeextrap_max = max(max(this.s3_nodeextrap));
            end
            if (this.taumax)
                this.tmax_gp_min         = min(min(this.tmax_gp));
                this.tmax_gp_max         = max(max(this.tmax_gp));
                this.tmax_elemextrap_min = min(min(this.tmax_elemextrap));
                this.tmax_elemextrap_max = max(max(this.tmax_elemextrap));
                this.tmax_nodeextrap_min = min(min(this.tmax_nodeextrap));
                this.tmax_nodeextrap_max = max(max(this.tmax_nodeextrap));
            end
            if (this.fxx)
                this.fxx_gp_min         = min(min(this.fxx_gp));
                this.fxx_gp_max         = max(max(this.fxx_gp));
                this.fxx_elemextrap_min = min(min(this.fxx_elemextrap));
                this.fxx_elemextrap_max = max(max(this.fxx_elemextrap));
                this.fxx_nodeextrap_min = min(min(this.fxx_nodeextrap));
                this.fxx_nodeextrap_max = max(max(this.fxx_nodeextrap));
            end
            if (this.fyy)
                this.fyy_gp_min         = min(min(this.fyy_gp));
                this.fyy_gp_max         = max(max(this.fyy_gp));
                this.fyy_elemextrap_min = min(min(this.fyy_elemextrap));
                this.fyy_elemextrap_max = max(max(this.fyy_elemextrap));
                this.fyy_nodeextrap_min = min(min(this.fyy_nodeextrap));
                this.fyy_nodeextrap_max = max(max(this.fyy_nodeextrap));
            end
            if (this.fzz)
                this.fzz_gp_min         = min(min(this.fzz_gp));
                this.fzz_gp_max         = max(max(this.fzz_gp));
                this.fzz_elemextrap_min = min(min(this.fzz_elemextrap));
                this.fzz_elemextrap_max = max(max(this.fzz_elemextrap));
                this.fzz_nodeextrap_min = min(min(this.fzz_nodeextrap));
                this.fzz_nodeextrap_max = max(max(this.fzz_nodeextrap));
            end
            if (this.fm)
                this.fm_gp_min         = min(min(this.fm_gp));
                this.fm_gp_max         = max(max(this.fm_gp));
                this.fm_elemextrap_min = min(min(this.fm_elemextrap));
                this.fm_elemextrap_max = max(max(this.fm_elemextrap));
                this.fm_nodeextrap_min = min(min(this.fm_nodeextrap));
                this.fm_nodeextrap_max = max(max(this.fm_nodeextrap));
            end
        end
        
        %------------------------------------------------------------------
        % Clear numerical garbage from pos-process results and smooth
        % constant results. (TEMPORARY SIMPLIFICATION !!)
        function clearSmallValues(this,mdl)
            if (this.sxx)
                % Clear small values
                if (abs(this.sxx_gp_min) < this.tol && abs(this.sxx_gp_max) < this.tol)
                    this.sxx_gp_min = 0.0;
                    this.sxx_gp_max = 0.0;
                    this.sxx_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.sxx_elemextrap_min) < this.tol && abs(this.sxx_elemextrap_max) < this.tol)
                    this.sxx_elemextrap_min = 0.0;
                    this.sxx_elemextrap_max = 0.0;
                    this.sxx_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.sxx_nodeextrap_min) < this.tol && abs(this.sxx_nodeextrap_max) < this.tol)
                    this.sxx_nodeextrap_min = 0.0;
                    this.sxx_nodeextrap_max = 0.0;
                    this.sxx_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.sxx_gp_min - this.sxx_gp_max) < this.tol)
                    mean = (this.sxx_gp_min+this.sxx_gp_max)/2;
                    this.sxx_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.sxx_gp_min = this.sxx_gp_min-abs(this.sxx_gp_min/10);
                    this.sxx_gp_max = this.sxx_gp_max+abs(this.sxx_gp_max/10);
                end
                if (abs(this.sxx_elemextrap_min - this.sxx_elemextrap_max) < this.tol)
                    mean = (this.sxx_elemextrap_min+this.sxx_elemextrap_max)/2;
                    this.sxx_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.sxx_elemextrap_min = this.sxx_elemextrap_min-abs(this.sxx_elemextrap_min/10);
                    this.sxx_elemextrap_max = this.sxx_elemextrap_max+abs(this.sxx_elemextrap_max/10);
                end
                if (abs(this.sxx_nodeextrap_min - this.sxx_nodeextrap_max) < this.tol)
                    mean = (this.sxx_nodeextrap_min+this.sxx_nodeextrap_max)/2;
                    this.sxx_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.sxx_nodeextrap_min = this.sxx_nodeextrap_min-abs(this.sxx_nodeextrap_min/10);
                    this.sxx_nodeextrap_max = this.sxx_nodeextrap_max+abs(this.sxx_nodeextrap_max/10);
                end
            end
            if (this.syy)
                % Clear small values
                if (abs(this.syy_gp_min) < this.tol && abs(this.syy_gp_max) < this.tol)
                    this.syy_gp_min = 0.0;
                    this.syy_gp_max = 0.0;
                    this.syy_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.syy_elemextrap_min) < this.tol && abs(this.syy_elemextrap_max) < this.tol)
                    this.syy_elemextrap_min = 0.0;
                    this.syy_elemextrap_max = 0.0;
                    this.syy_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.syy_nodeextrap_min) < this.tol && abs(this.syy_nodeextrap_max) < this.tol)
                    this.syy_nodeextrap_min = 0.0;
                    this.syy_nodeextrap_max = 0.0;
                    this.syy_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.syy_gp_min - this.syy_gp_max) < this.tol)
                    mean = (this.syy_gp_min+this.syy_gp_max)/2;
                    this.syy_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.syy_gp_min = this.syy_gp_min-abs(this.syy_gp_min/10);
                    this.syy_gp_max = this.syy_gp_max+abs(this.syy_gp_max/10);
                end
                if (abs(this.syy_elemextrap_min - this.syy_elemextrap_max) < this.tol)
                    mean = (this.syy_elemextrap_min+this.syy_elemextrap_max)/2;
                    this.syy_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.syy_elemextrap_min = this.syy_elemextrap_min-abs(this.syy_elemextrap_min/10);
                    this.syy_elemextrap_max = this.syy_elemextrap_max+abs(this.syy_elemextrap_max/10);
                end
                if (abs(this.syy_nodeextrap_min - this.syy_nodeextrap_max) < this.tol)
                    mean = (this.syy_nodeextrap_min+this.syy_nodeextrap_max)/2;
                    this.syy_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.syy_nodeextrap_min = this.syy_nodeextrap_min-abs(this.syy_nodeextrap_min/10);
                    this.syy_nodeextrap_max = this.syy_nodeextrap_max+abs(this.syy_nodeextrap_max/10);
                end
            end
            if (this.szz)
                % Clear small values
                if (abs(this.szz_gp_min) < this.tol && abs(this.szz_gp_max) < this.tol)
                    this.szz_gp_min = 0.0;
                    this.szz_gp_max = 0.0;
                    this.szz_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.szz_elemextrap_min) < this.tol && abs(this.szz_elemextrap_max) < this.tol)
                    this.szz_elemextrap_min = 0.0;
                    this.szz_elemextrap_max = 0.0;
                    this.szz_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.szz_nodeextrap_min) < this.tol && abs(this.szz_nodeextrap_max) < this.tol)
                    this.szz_nodeextrap_min = 0.0;
                    this.szz_nodeextrap_max = 0.0;
                    this.szz_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.szz_gp_min - this.szz_gp_max) < this.tol)
                    mean = (this.szz_gp_min+this.szz_gp_max)/2;
                    this.szz_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.szz_gp_min = this.szz_gp_min-abs(this.szz_gp_min/10);
                    this.szz_gp_max = this.szz_gp_max+abs(this.szz_gp_max/10);
                end
                if (abs(this.szz_elemextrap_min - this.szz_elemextrap_max) < this.tol)
                    mean = (this.szz_elemextrap_min+this.szz_elemextrap_max)/2;
                    this.szz_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.szz_elemextrap_min = this.szz_elemextrap_min-abs(this.szz_elemextrap_min/10);
                    this.szz_elemextrap_max = this.szz_elemextrap_max+abs(this.szz_elemextrap_max/10);
                end
                if (abs(this.szz_nodeextrap_min - this.szz_nodeextrap_max) < this.tol)
                    mean = (this.szz_nodeextrap_min+this.szz_nodeextrap_max)/2;
                    this.szz_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.szz_nodeextrap_min = this.szz_nodeextrap_min-abs(this.szz_nodeextrap_min/10);
                    this.szz_nodeextrap_max = this.szz_nodeextrap_max+abs(this.szz_nodeextrap_max/10);
                end
            end
            if (this.txy)
                % Clear small values
                if (abs(this.txy_gp_min) < this.tol && abs(this.txy_gp_max) < this.tol)
                    this.txy_gp_min = 0.0;
                    this.txy_gp_max = 0.0;
                    this.txy_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.txy_elemextrap_min) < this.tol && abs(this.txy_elemextrap_max) < this.tol)
                    this.txy_elemextrap_min = 0.0;
                    this.txy_elemextrap_max = 0.0;
                    this.txy_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.txy_nodeextrap_min) < this.tol && abs(this.txy_nodeextrap_max) < this.tol)
                    this.txy_nodeextrap_min = 0.0;
                    this.txy_nodeextrap_max = 0.0;
                    this.txy_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.txy_gp_min - this.txy_gp_max) < this.tol)
                    mean = (this.txy_gp_min+this.txy_gp_max)/2;
                    this.txy_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.txy_gp_min = this.txy_gp_min-abs(this.txy_gp_min/10);
                    this.txy_gp_max = this.txy_gp_max+abs(this.txy_gp_max/10);
                end
                if (abs(this.txy_elemextrap_min - this.txy_elemextrap_max) < this.tol)
                    mean = (this.txy_elemextrap_min+this.txy_elemextrap_max)/2;
                    this.txy_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.txy_elemextrap_min = this.txy_elemextrap_min-abs(this.txy_elemextrap_min/10);
                    this.txy_elemextrap_max = this.txy_elemextrap_max+abs(this.txy_elemextrap_max/10);
                end
                if (abs(this.txy_nodeextrap_min - this.txy_nodeextrap_max) < this.tol)
                    mean = (this.txy_nodeextrap_min+this.txy_nodeextrap_max)/2;
                    this.txy_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.txy_nodeextrap_min = this.txy_nodeextrap_min-abs(this.txy_nodeextrap_min/10);
                    this.txy_nodeextrap_max = this.txy_nodeextrap_max+abs(this.txy_nodeextrap_max/10);
                end
            end
            if (this.txz)
                % Clear small values
                if (abs(this.txz_gp_min) < this.tol && abs(this.txz_gp_max) < this.tol)
                    this.txz_gp_min = 0.0;
                    this.txz_gp_max = 0.0;
                    this.txz_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.txz_elemextrap_min) < this.tol && abs(this.txz_elemextrap_max) < this.tol)
                    this.txz_elemextrap_min = 0.0;
                    this.txz_elemextrap_max = 0.0;
                    this.txz_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.txz_nodeextrap_min) < this.tol && abs(this.txz_nodeextrap_max) < this.tol)
                    this.txz_nodeextrap_min = 0.0;
                    this.txz_nodeextrap_max = 0.0;
                    this.txz_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.txz_gp_min - this.txz_gp_max) < this.tol)
                    mean = (this.txz_gp_min+this.txz_gp_max)/2;
                    this.txz_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.txz_gp_min = this.txz_gp_min-abs(this.txz_gp_min/10);
                    this.txz_gp_max = this.txz_gp_max+abs(this.txz_gp_max/10);
                end
                if (abs(this.txz_elemextrap_min - this.txz_elemextrap_max) < this.tol)
                    mean = (this.txz_elemextrap_min+this.txz_elemextrap_max)/2;
                    this.txz_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.txz_elemextrap_min = this.txz_elemextrap_min-abs(this.txz_elemextrap_min/10);
                    this.txz_elemextrap_max = this.txz_elemextrap_max+abs(this.txz_elemextrap_max/10);
                end
                if (abs(this.txz_nodeextrap_min - this.txz_nodeextrap_max) < this.tol)
                    mean = (this.txz_nodeextrap_min+this.txz_nodeextrap_max)/2;
                    this.txz_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.txz_nodeextrap_min = this.txz_nodeextrap_min-abs(this.txz_nodeextrap_min/10);
                    this.txz_nodeextrap_max = this.txz_nodeextrap_max+abs(this.txz_nodeextrap_max/10);
                end
            end
            if (this.tyz)
                % Clear small values
                if (abs(this.tyz_gp_min) < this.tol && abs(this.tyz_gp_max) < this.tol)
                    this.tyz_gp_min = 0.0;
                    this.tyz_gp_max = 0.0;
                    this.tyz_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.tyz_elemextrap_min) < this.tol && abs(this.tyz_elemextrap_max) < this.tol)
                    this.tyz_elemextrap_min = 0.0;
                    this.tyz_elemextrap_max = 0.0;
                    this.tyz_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.tyz_nodeextrap_min) < this.tol && abs(this.tyz_nodeextrap_max) < this.tol)
                    this.tyz_nodeextrap_min = 0.0;
                    this.tyz_nodeextrap_max = 0.0;
                    this.tyz_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.tyz_gp_min - this.tyz_gp_max) < this.tol)
                    mean = (this.tyz_gp_min+this.tyz_gp_max)/2;
                    this.tyz_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.tyz_gp_min = this.tyz_gp_min-abs(this.tyz_gp_min/10);
                    this.tyz_gp_max = this.tyz_gp_max+abs(this.tyz_gp_max/10);
                end
                if (abs(this.tyz_elemextrap_min - this.tyz_elemextrap_max) < this.tol)
                    mean = (this.tyz_elemextrap_min+this.tyz_elemextrap_max)/2;
                    this.tyz_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.tyz_elemextrap_min = this.tyz_elemextrap_min-abs(this.tyz_elemextrap_min/10);
                    this.tyz_elemextrap_max = this.tyz_elemextrap_max+abs(this.tyz_elemextrap_max/10);
                end
                if (abs(this.tyz_nodeextrap_min - this.tyz_nodeextrap_max) < this.tol)
                    mean = (this.tyz_nodeextrap_min+this.tyz_nodeextrap_max)/2;
                    this.tyz_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.tyz_nodeextrap_min = this.tyz_nodeextrap_min-abs(this.tyz_nodeextrap_min/10);
                    this.tyz_nodeextrap_max = this.tyz_nodeextrap_max+abs(this.tyz_nodeextrap_max/10);
                end
            end
            if (this.s1)
                % Clear small values
                if (abs(this.s1_gp_min) < this.tol && abs(this.s1_gp_max) < this.tol)
                    this.s1_gp_min = 0.0;
                    this.s1_gp_max = 0.0;
                    this.s1_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.s1_elemextrap_min) < this.tol && abs(this.s1_elemextrap_max) < this.tol)
                    this.s1_elemextrap_min = 0.0;
                    this.s1_elemextrap_max = 0.0;
                    this.s1_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.s1_nodeextrap_min) < this.tol && abs(this.s1_nodeextrap_max) < this.tol)
                    this.s1_nodeextrap_min = 0.0;
                    this.s1_nodeextrap_max = 0.0;
                    this.s1_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.s1_gp_min - this.s1_gp_max) < this.tol)
                    mean = (this.s1_gp_min+this.s1_gp_max)/2;
                    this.s1_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.s1_gp_min = this.s1_gp_min-abs(this.s1_gp_min/10);
                    this.s1_gp_max = this.s1_gp_max+abs(this.s1_gp_max/10);
                end
                if (abs(this.s1_elemextrap_min - this.s1_elemextrap_max) < this.tol)
                    mean = (this.s1_elemextrap_min+this.s1_elemextrap_max)/2;
                    this.s1_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.s1_elemextrap_min = this.s1_elemextrap_min-abs(this.s1_elemextrap_min/10);
                    this.s1_elemextrap_max = this.s1_elemextrap_max+abs(this.s1_elemextrap_max/10);
                end
                if (abs(this.s1_nodeextrap_min - this.s1_nodeextrap_max) < this.tol)
                    mean = (this.s1_nodeextrap_min+this.s1_nodeextrap_max)/2;
                    this.s1_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.s1_nodeextrap_min = this.s1_nodeextrap_min-abs(this.s1_nodeextrap_min/10);
                    this.s1_nodeextrap_max = this.s1_nodeextrap_max+abs(this.s1_nodeextrap_max/10);
                end
            end
            if (this.s2)
                % Clear small values
                if (abs(this.s2_gp_min) < this.tol && abs(this.s2_gp_max) < this.tol)
                    this.s2_gp_min = 0.0;
                    this.s2_gp_max = 0.0;
                    this.s2_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.s2_elemextrap_min) < this.tol && abs(this.s2_elemextrap_max) < this.tol)
                    this.s2_elemextrap_min = 0.0;
                    this.s2_elemextrap_max = 0.0;
                    this.s2_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.s2_nodeextrap_min) < this.tol && abs(this.s2_nodeextrap_max) < this.tol)
                    this.s2_nodeextrap_min = 0.0;
                    this.s2_nodeextrap_max = 0.0;
                    this.s2_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.s2_gp_min - this.s2_gp_max) < this.tol)
                    mean = (this.s2_gp_min+this.s2_gp_max)/2;
                    this.s2_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.s2_gp_min = this.s2_gp_min-abs(this.s2_gp_min/10);
                    this.s2_gp_max = this.s2_gp_max+abs(this.s2_gp_max/10);
                end
                if (abs(this.s2_elemextrap_min - this.s2_elemextrap_max) < this.tol)
                    mean = (this.s2_elemextrap_min+this.s2_elemextrap_max)/2;
                    this.s2_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.s2_elemextrap_min = this.s2_elemextrap_min-abs(this.s2_elemextrap_min/10);
                    this.s2_elemextrap_max = this.s2_elemextrap_max+abs(this.s2_elemextrap_max/10);
                end
                if (abs(this.s2_nodeextrap_min - this.s2_nodeextrap_max) < this.tol)
                    mean = (this.s2_nodeextrap_min+this.s2_nodeextrap_max)/2;
                    this.s2_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.s2_nodeextrap_min = this.s2_nodeextrap_min-abs(this.s2_nodeextrap_min/10);
                    this.s2_nodeextrap_max = this.s2_nodeextrap_max+abs(this.s2_nodeextrap_max/10);
                end
            end
            if (this.s3)
                % Clear small values
                if (abs(this.s3_gp_min) < this.tol && abs(this.s3_gp_max) < this.tol)
                    this.s3_gp_min = 0.0;
                    this.s3_gp_max = 0.0;
                    this.s3_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.s3_elemextrap_min) < this.tol && abs(this.S3_elemextrap_max) < this.tol)
                    this.s3_elemextrap_min = 0.0;
                    this.s3_elemextrap_max = 0.0;
                    this.s3_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.s3_nodeextrap_min) < this.tol && abs(this.s3_nodeextrap_max) < this.tol)
                    this.s3_nodeextrap_min = 0.0;
                    this.s3_nodeextrap_max = 0.0;
                    this.s3_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.s3_gp_min - this.s3_gp_max) < this.tol)
                    mean = (this.s3_gp_min+this.s3_gp_max)/2;
                    this.s3_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.s3_gp_min = this.s3_gp_min-abs(this.s3_gp_min/10);
                    this.s3_gp_max = this.s3_gp_max+abs(this.s3_gp_max/10);
                end
                if (abs(this.s3_elemextrap_min - this.s3_elemextrap_max) < this.tol)
                    mean = (this.s3_elemextrap_min+this.s3_elemextrap_max)/2;
                    this.s3_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.s3_elemextrap_min = this.s3_elemextrap_min-abs(this.s3_elemextrap_min/10);
                    this.s3_elemextrap_max = this.s3_elemextrap_max+abs(this.s3_elemextrap_max/10);
                end
                if (abs(this.s3_nodeextrap_min - this.s3_nodeextrap_max) < this.tol)
                    mean = (this.s3_nodeextrap_min+this.s3_nodeextrap_max)/2;
                    this.s3_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.s3_nodeextrap_min = this.s3_nodeextrap_min-abs(this.s3_nodeextrap_min/10);
                    this.s3_nodeextrap_max = this.s3_nodeextrap_max+abs(this.s3_nodeextrap_max/10);
                end
            end
            if (this.taumax)
                % Clear small values
                if (abs(this.tmax_gp_min) < this.tol && abs(this.tmax_gp_max) < this.tol)
                    this.tmax_gp_min = 0.0;
                    this.tmax_gp_max = 0.0;
                    this.tmax_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.tmax_elemextrap_min) < this.tol && abs(this.tmax_elemextrap_max) < this.tol)
                    this.tmax_elemextrap_min = 0.0;
                    this.tmax_elemextrap_max = 0.0;
                    this.tmax_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.tmax_nodeextrap_min) < this.tol && abs(this.tmax_nodeextrap_max) < this.tol)
                    this.tmax_nodeextrap_min = 0.0;
                    this.tmax_nodeextrap_max = 0.0;
                    this.tmax_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.tmax_gp_min - this.tmax_gp_max) < this.tol)
                    mean = (this.tmax_gp_min+this.tmax_gp_max)/2;
                    this.tmax_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.tmax_gp_min = this.tmax_gp_min-abs(this.tmax_gp_min/10);
                    this.tmax_gp_max = this.tmax_gp_max+abs(this.tmax_gp_max/10);
                end
                if (abs(this.tmax_elemextrap_min - this.tmax_elemextrap_max) < this.tol)
                    mean = (this.tmax_elemextrap_min+this.tmax_elemextrap_max)/2;
                    this.tmax_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.tmax_elemextrap_min = this.tmax_elemextrap_min-abs(this.tmax_elemextrap_min/10);
                    this.tmax_elemextrap_max = this.tmax_elemextrap_max+abs(this.tmax_elemextrap_max/10);
                end
                if (abs(this.tmax_nodeextrap_min - this.tmax_nodeextrap_max) < this.tol)
                    mean = (this.tmax_nodeextrap_min+this.tmax_nodeextrap_max)/2;
                    this.tmax_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.tmax_nodeextrap_min = this.tmax_nodeextrap_min-abs(this.tmax_nodeextrap_min/10);
                    this.tmax_nodeextrap_max = this.tmax_nodeextrap_max+abs(this.tmax_nodeextrap_max/10);
                end
            end
            if (this.fxx)
                % Clear small values
                if (abs(this.fxx_gp_min) < this.tol && abs(this.fxx_gp_max) < this.tol)
                    this.fxx_gp_min = 0.0;
                    this.fxx_gp_max = 0.0;
                    this.fxx_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.fxx_elemextrap_min) < this.tol && abs(this.fxx_elemextrap_max) < this.tol)
                    this.fxx_elemextrap_min = 0.0;
                    this.fxx_elemextrap_max = 0.0;
                    this.fxx_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.fxx_nodeextrap_min) < this.tol && abs(this.fxx_nodeextrap_max) < this.tol)
                    this.fxx_nodeextrap_min = 0.0;
                    this.fxx_nodeextrap_max = 0.0;
                    this.fxx_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.fxx_gp_min - this.fxx_gp_max) < this.tol)
                    mean = (this.fxx_gp_min+this.fxx_gp_max)/2;
                    this.fxx_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.fxx_gp_min = this.fxx_gp_min-abs(this.fxx_gp_min/10);
                    this.fxx_gp_max = this.fxx_gp_max+abs(this.fxx_gp_max/10);
                end
                if (abs(this.fxx_elemextrap_min - this.fxx_elemextrap_max) < this.tol)
                    mean = (this.fxx_elemextrap_min+this.fxx_elemextrap_max)/2;
                    this.fxx_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.fxx_elemextrap_min = this.fxx_elemextrap_min-abs(this.fxx_elemextrap_min/10);
                    this.fxx_elemextrap_max = this.fxx_elemextrap_max+abs(this.fxx_elemextrap_max/10);
                end
                if (abs(this.fxx_nodeextrap_min - this.fxx_nodeextrap_max) < this.tol)
                    mean = (this.fxx_nodeextrap_min+this.fxx_nodeextrap_max)/2;
                    this.fxx_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.fxx_nodeextrap_min = this.fxx_nodeextrap_min-abs(this.fxx_nodeextrap_min/10);
                    this.fxx_nodeextrap_max = this.fxx_nodeextrap_max+abs(this.fxx_nodeextrap_max/10);
                end
            end
            if (this.fyy)
                % Clear small values
                if (abs(this.fyy_gp_min) < this.tol && abs(this.fyy_gp_max) < this.tol)
                    this.fyy_gp_min = 0.0;
                    this.fyy_gp_max = 0.0;
                    this.fyy_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.fyy_elemextrap_min) < this.tol && abs(this.fyy_elemextrap_max) < this.tol)
                    this.fyy_elemextrap_min = 0.0;
                    this.fyy_elemextrap_max = 0.0;
                    this.fyy_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.fyy_nodeextrap_min) < this.tol && abs(this.fyy_nodeextrap_max) < this.tol)
                    this.fyy_nodeextrap_min = 0.0;
                    this.fyy_nodeextrap_max = 0.0;
                    this.fyy_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.fyy_gp_min - this.fyy_gp_max) < this.tol)
                    mean = (this.fyy_gp_min+this.fyy_gp_max)/2;
                    this.fyy_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.fyy_gp_min = this.fyy_gp_min-abs(this.fyy_gp_min/10);
                    this.fyy_gp_max = this.fyy_gp_max+abs(this.fyy_gp_max/10);
                end
                if (abs(this.fyy_elemextrap_min - this.fyy_elemextrap_max) < this.tol)
                    mean = (this.fyy_elemextrap_min+this.fyy_elemextrap_max)/2;
                    this.fyy_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.fyy_elemextrap_min = this.fyy_elemextrap_min-abs(this.fyy_elemextrap_min/10);
                    this.fyy_elemextrap_max = this.fyy_elemextrap_max+abs(this.fyy_elemextrap_max/10);
                end
                if (abs(this.fyy_nodeextrap_min - this.fyy_nodeextrap_max) < this.tol)
                    mean = (this.fyy_nodeextrap_min+this.fyy_nodeextrap_max)/2;
                    this.fyy_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.fyy_nodeextrap_min = this.fyy_nodeextrap_min-abs(this.fyy_nodeextrap_min/10);
                    this.fyy_nodeextrap_max = this.fyy_nodeextrap_max+abs(this.fyy_nodeextrap_max/10);
                end
            end
            if (this.fzz)
                % Clear small values
                if (abs(this.fzz_gp_min) < this.tol && abs(this.fzz_gp_max) < this.tol)
                    this.fzz_gp_min = 0.0;
                    this.fzz_gp_max = 0.0;
                    this.fzz_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.fzz_elemextrap_min) < this.tol && abs(this.fzz_elemextrap_max) < this.tol)
                    this.fzz_elemextrap_min = 0.0;
                    this.fzz_elemextrap_max = 0.0;
                    this.fzz_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.fzz_nodeextrap_min) < this.tol && abs(this.fzz_nodeextrap_max) < this.tol)
                    this.fzz_nodeextrap_min = 0.0;
                    this.fzz_nodeextrap_max = 0.0;
                    this.fzz_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.fzz_gp_min - this.fzz_gp_max) < this.tol)
                    mean = (this.fzz_gp_min+this.fzz_gp_max)/2;
                    this.fzz_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.fzz_gp_min = this.fzz_gp_min-abs(this.fzz_gp_min/10);
                    this.fzz_gp_max = this.fzz_gp_max+abs(this.fzz_gp_max/10);
                end
                if (abs(this.fzz_elemextrap_min - this.fzz_elemextrap_max) < this.tol)
                    mean = (this.fzz_elemextrap_min+this.fzz_elemextrap_max)/2;
                    this.fzz_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.fzz_elemextrap_min = this.fzz_elemextrap_min-abs(this.fzz_elemextrap_min/10);
                    this.fzz_elemextrap_max = this.fzz_elemextrap_max+abs(this.fzz_elemextrap_max/10);
                end
                if (abs(this.fzz_nodeextrap_min - this.fzz_nodeextrap_max) < this.tol)
                    mean = (this.fzz_nodeextrap_min+this.fzz_nodeextrap_max)/2;
                    this.fzz_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.fzz_nodeextrap_min = this.fzz_nodeextrap_min-abs(this.fzz_nodeextrap_min/10);
                    this.fzz_nodeextrap_max = this.fzz_nodeextrap_max+abs(this.fzz_nodeextrap_max/10);
                end
            end
            if (this.fm)
                % Clear small values
                if (abs(this.fm_gp_min) < this.tol && abs(this.fm_gp_max) < this.tol)
                    this.fm_gp_min = 0.0;
                    this.fm_gp_max = 0.0;
                    this.fm_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.fm_elemextrap_min) < this.tol && abs(this.fm_elemextrap_max) < this.tol)
                    this.fm_elemextrap_min = 0.0;
                    this.fm_elemextrap_max = 0.0;
                    this.fm_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.fm_nodeextrap_min) < this.tol && abs(this.fm_nodeextrap_max) < this.tol)
                    this.fm_nodeextrap_min = 0.0;
                    this.fm_nodeextrap_max = 0.0;
                    this.fm_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.fm_gp_min - this.fm_gp_max) < this.tol)
                    mean = (this.fm_gp_min+this.fm_gp_max)/2;
                    this.fm_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.fm_gp_min = this.fm_gp_min-abs(this.fm_gp_min/10);
                    this.fm_gp_max = this.fm_gp_max+abs(this.fm_gp_max/10);
                end
                if (abs(this.fm_elemextrap_min - this.fm_elemextrap_max) < this.tol)
                    mean = (this.fm_elemextrap_min+this.fm_elemextrap_max)/2;
                    this.fm_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.fm_elemextrap_min = this.fm_elemextrap_min-abs(this.fm_elemextrap_min/10);
                    this.fm_elemextrap_max = this.fm_elemextrap_max+abs(this.fm_elemextrap_max/10);
                end
                if (abs(this.fm_nodeextrap_min - this.fm_nodeextrap_max) < this.tol)
                    mean = (this.fm_nodeextrap_min+this.fm_nodeextrap_max)/2;
                    this.fm_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.fm_nodeextrap_min = this.fm_nodeextrap_min-abs(this.fm_nodeextrap_min/10);
                    this.fm_nodeextrap_max = this.fm_nodeextrap_max+abs(this.fm_nodeextrap_max/10);
                end
            end
        end
    end
end
