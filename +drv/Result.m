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
        qxz    = logical(false);  % XZ shear force
        qyz    = logical(false);  % YZ shear force
        mxx    = logical(false);  % moment about X direction
        myy    = logical(false);  % moment about Y direction
        mxy    = logical(false);  % torsion moment
        m1     = logical(false);  % principal moment 1
        m2     = logical(false);  % principal moment 2
        tormax = logical(false);  % maximum torsion
        
        % Thermal analysis contours
        temp   = logical(false);  % temperature field
        fxx    = logical(false);  % heat fluxes in X direction
        fyy    = logical(false);  % heat fluxes in Y direction
        fzz    = logical(false);  % heat fluxes in Z direction
        fm     = logical(false);  % heat fluxes module
        
        % Result curves options
        curve_temp  int32 = int32.empty; % vector of node IDs to plot transient temperature response
        output_freq = 1;                 % output frequency for plotting transient response
        
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Results from equilibrium system
        U     double = double.empty;  % global vector of state variables
        Ut    double = double.empty;  % global vector of state variable 1st time derivatives
        Utt   double = double.empty;  % global vector of state variable 2nd time derivatives
        steps int32  = int32.empty;   % number of performed steps
        times double = double.empty;  % vector of time values of each step
        
        % Reference values to dimension pos-process arrays
        maxGPts int32 = int32.empty;  % max. number of Gauss points in all elements
        maxNen  int32 = int32.empty;  % max. number of nodal points in all elements
        
        % General Gauss point data
        ngp  int32  = int32.empty;  % vector of number of element gauss pts
        x_gp double = double.empty; % vector of gauss point x coordinates
        y_gp double = double.empty; % vector of gauss point y coordinates
        z_gp double = double.empty; % vector of gauss point Z coordinates
        
        % Structural analysis stresses at Gauss points
        sxx_gp                double = double.empty; % sigma x gauss points stress array
        sxx_gp_min            double = double.empty; % sigma x gauss points min. stress
        sxx_gp_max            double = double.empty; % sigma x gauss points max. stress
        syy_gp                double = double.empty; % sigma y gauss points stress array
        syy_gp_min            double = double.empty; % sigma y gauss points min. stress
        syy_gp_max            double = double.empty; % sigma y gauss points max. stress
        szz_gp                double = double.empty; % sigma Z gauss points stress array
        szz_gp_min            double = double.empty; % sigma Z gauss points min. stress
        szz_gp_max            double = double.empty; % sigma Z gauss points max. stress
        txy_gp                double = double.empty; % tau xy gauss points stress array
        txy_gp_min            double = double.empty; % tau xy gauss points min. stress
        txy_gp_max            double = double.empty; % tau xy gauss points max. stress
        txz_gp                double = double.empty; % tau xz gauss points stress array
        txz_gp_min            double = double.empty; % tau xz gauss points min. stress
        txz_gp_max            double = double.empty; % tau xz gauss points max. stress
        tyz_gp                double = double.empty; % tau yz gauss points stress array
        tyz_gp_min            double = double.empty; % tau yz gauss points min. stress
        tyz_gp_max            double = double.empty; % tau yz gauss points max. stress
        s1_gp                 double = double.empty; % sigma 1 gauss points stress array
        s1_gp_min             double = double.empty; % sigma 1 gauss points min. stress
        s1_gp_max             double = double.empty; % sigma 1 gauss points max. stress
        s2_gp                 double = double.empty; % sigma 2 gauss points stress array
        s2_gp_min             double = double.empty; % sigma 2 gauss points min. stress
        s2_gp_max             double = double.empty; % sigma 2 gauss points max. stress
        s3_gp                 double = double.empty; % sigma 3 gauss points stress array
        s3_gp_min             double = double.empty; % sigma 3 gauss points min. stress
        s3_gp_max             double = double.empty; % sigma 3 gauss points max. stress
        tmax_gp               double = double.empty; % tau max. gauss points stress array
        tmax_gp_min           double = double.empty; % tau max. gauss points min. stress
        tmax_gp_max           double = double.empty; % tau max. gauss points max. stress
        s1x_gp                double = double.empty; % vector of s1 vector gauss x components
        s1y_gp                double = double.empty; % vector of s1 vector gauss y components
        s1z_gp                double = double.empty; % vector of s1 vector gauss z components
        s2x_gp                double = double.empty; % vector of s2 vector gauss x components
        s2y_gp                double = double.empty; % vector of s2 vector gauss y components
        s2z_gp                double = double.empty; % vector of s2 vector gauss y components
        s3x_gp                double = double.empty; % vector of s3 vector gauss x components
        s3y_gp                double = double.empty; % vector of s3 vector gauss y components
        s3z_gp                double = double.empty; % vector of s3 vector gauss z components
        
        % Structural analysis internal forces at Gauss points
        qxz_gp                double = double.empty; % shear xz gauss points array
        qxz_gp_min            double = double.empty; % shear xz gauss points min. value
        qxz_gp_max            double = double.empty; % shear xz gauss points max. value
        qyz_gp                double = double.empty; % shear yz gauss points array
        qyz_gp_min            double = double.empty; % shear yz gauss points min. value
        qyz_gp_max            double = double.empty; % shear yz gauss points max. value
        mxx_gp                double = double.empty; % moment x gauss points array
        mxx_gp_min            double = double.empty; % moment x gauss points min. value
        mxx_gp_max            double = double.empty; % moment x gauss points max. value
        myy_gp                double = double.empty; % moment y gauss points array
        myy_gp_min            double = double.empty; % moment y gauss points min. value
        myy_gp_max            double = double.empty; % moment y gauss points max. value
        mxy_gp                double = double.empty; % moment xy gauss points array
        mxy_gp_min            double = double.empty; % moment xy gauss points min. value
        mxy_gp_max            double = double.empty; % moment xy gauss points max. value
        m1_gp                 double = double.empty; % moment 1 gauss points array
        m1_gp_min             double = double.empty; % moment 1 gauss points min. value
        m1_gp_max             double = double.empty; % moment 1 gauss points max. value
        m2_gp                 double = double.empty; % moment 2 gauss points array
        m2_gp_min             double = double.empty; % moment 2 gauss points min. value
        m2_gp_max             double = double.empty; % moment 2 gauss points max. value
        m1x_gp                double = double.empty; % vector of gauss moment 1 x components
        m1y_gp                double = double.empty; % vector of gauss moment 1 y components
        m2x_gp                double = double.empty; % vector of gauss moment 2 x components
        m2y_gp                double = double.empty; % vector of gauss moment 2 y components
        tormax_gp             double = double.empty; % torsion max. gauss points array
        tormax_gp_min         double = double.empty; % torsion max. gauss points min. value
        tormax_gp_max         double = double.empty; % torsion max. gauss points max. value
        
        % Thermal analysis heat fluxes at Gauss points
        fxx_gp                double = double.empty; % heat flux x gauss points array
        fxx_gp_min            double = double.empty; % heat flux x gauss points min. value
        fxx_gp_max            double = double.empty; % heat flux x gauss points max. value
        fyy_gp                double = double.empty; % heat flux y gauss points array
        fyy_gp_min            double = double.empty; % heat flux y gauss points min. value
        fyy_gp_max            double = double.empty; % heat flux y gauss points max. value
        fzz_gp                double = double.empty; % heat flux z gauss points array
        fzz_gp_min            double = double.empty; % heat flux z gauss points min. value
        fzz_gp_max            double = double.empty; % heat flux z gauss points max. value
        fm_gp                 double = double.empty; % heat flux module gauss points array
        fm_gp_min             double = double.empty; % heat flux module gauss points min. value
        fm_gp_max             double = double.empty; % heat flux module gauss points max. value
        fmx_gp                double = double.empty; % vector of gauss heat flux module x components
        fmy_gp                double = double.empty; % vector of gauss heat flux module y components
        fmz_gp                double = double.empty; % vector of gauss heat flux module z components
        
        % Structural analysis stresses extrapolated to element nodes
        sxx_elemextrap        double = double.empty; % sigma x element node extrap. stress array
        sxx_elemextrap_min    double = double.empty; % sigma x element node extrap. min. stress
        sxx_elemextrap_max    double = double.empty; % sigma x element node extrap. max. stress
        syy_elemextrap        double = double.empty; % sigma y element node extrap. stress array
        syy_elemextrap_min    double = double.empty; % sigma y element node extrap. min. stress
        syy_elemextrap_max    double = double.empty; % sigma y element node extrap. max. stress
        szz_elemextrap        double = double.empty; % sigma z element node extrap. stress array
        szz_elemextrap_min    double = double.empty; % sigma z element node extrap. min. stress
        szz_elemextrap_max    double = double.empty; % sigma z element node extrap. max. stress
        txy_elemextrap        double = double.empty; % tau xy element node extrap. stress array
        txy_elemextrap_min    double = double.empty; % tau xy element node extrap. min. stress
        txy_elemextrap_max    double = double.empty; % tau xy element node extrap. max. stress
        txz_elemextrap        double = double.empty; % tau xz element node extrap. stress array
        txz_elemextrap_min    double = double.empty; % tau xz element node extrap. min. stress
        txz_elemextrap_max    double = double.empty; % tau xz element node extrap. max. stress
        tyz_elemextrap        double = double.empty; % tau yz element node extrap. stress array
        tyz_elemextrap_min    double = double.empty; % tau yz element node extrap. min. stress
        tyz_elemextrap_max    double = double.empty; % tau yz element node extrap. max. stress
        s1_elemextrap         double = double.empty; % sigma 1 element node extrap. stress array
        s1_elemextrap_min     double = double.empty; % sigma 1 element node extrap. min. stress
        s1_elemextrap_max     double = double.empty; % sigma 1 element node extrap. max. stress
        s2_elemextrap         double = double.empty; % sigma 2 element node extrap. stress array
        s2_elemextrap_min     double = double.empty; % sigma 2 element node extrap. min. stress
        s2_elemextrap_max     double = double.empty; % sigma 2 element node extrap. max. stress
        s3_elemextrap         double = double.empty; % sigma 3 element node extrap. stress array
        s3_elemextrap_min     double = double.empty; % sigma 3 element node extrap. min. stress
        s3_elemextrap_max     double = double.empty; % sigma 3 element node extrap. max. stress
        tmax_elemextrap       double = double.empty; % tau max. element node extrap. stress array
        tmax_elemextrap_min   double = double.empty; % tau max. element node extrap. min stress
        tmax_elemextrap_max   double = double.empty; % tau max. element node extrap. max stress
        
        % Structural analysis internal forces extrapolated to element nodes
        qxz_elemextrap        double = double.empty; % shear xz element node extrap. shear array
        qxz_elemextrap_min    double = double.empty; % shear xz element node extrap. min. shear array
        qxz_elemextrap_max    double = double.empty; % shear xz element node extrap. max. shear array
        qyz_elemextrap        double = double.empty; % shear yz element node extrap. shear array
        qyz_elemextrap_min    double = double.empty; % shear yz element node extrap. min. shear array
        qyz_elemextrap_max    double = double.empty; % shear yz element node extrap. max. shear array
        mxx_elemextrap        double = double.empty; % moment x element node extrap. moments array
        mxx_elemextrap_min    double = double.empty; % moment x element node extrap. min. moments array
        mxx_elemextrap_max    double = double.empty; % moment x element node extrap. max. moments array
        myy_elemextrap        double = double.empty; % moment y element node extrap. moments array
        myy_elemextrap_min    double = double.empty; % moment y element node extrap. min. moments array
        myy_elemextrap_max    double = double.empty; % moment y element node extrap. max. moments array
        mxy_elemextrap        double = double.empty; % moment xy element node extrap. moments array
        mxy_elemextrap_min    double = double.empty; % moment xy element node extrap. min. moments array
        mxy_elemextrap_max    double = double.empty; % moment xy element node extrap. max. moments array
        m1_elemextrap         double = double.empty; % moment 1 element node extrap. moments array
        m1_elemextrap_min     double = double.empty; % moment 1 element node extrap. min. moments array
        m1_elemextrap_max     double = double.empty; % moment 1 element node extrap. max. moments array
        m2_elemextrap         double = double.empty; % moment 2 element node extrap. moments array
        m2_elemextrap_min     double = double.empty; % moment 2 element node extrap. min. moments array
        m2_elemextrap_max     double = double.empty; % moment 2 element node extrap. max. moments array
        tormax_elemextrap     double = double.empty; % torsion max. element node extrap.  torsion array
        tormax_elemextrap_min double = double.empty; % torsion max. element node extrap. torsion array
        tormax_elemextrap_max double = double.empty; % torsion max. element node extrap. max. torsion array
        
        % Thermal analysis heat fluxes extrapolated to element nodes
        fxx_elemextrap        double = double.empty; % flux x element node extrap. array
        fxx_elemextrap_min    double = double.empty; % flux x element node extrap. min. value
        fxx_elemextrap_max    double = double.empty; % flux x element node extrap. max. value
        fyy_elemextrap        double = double.empty; % flux y element node extrap. array
        fyy_elemextrap_min    double = double.empty; % flux y element node extrap. min. value
        fyy_elemextrap_max    double = double.empty; % flux y element node extrap. max. value
        fzz_elemextrap        double = double.empty; % flux z element node extrap. array
        fzz_elemextrap_min    double = double.empty; % flux z element node extrap. min. value
        fzz_elemextrap_max    double = double.empty; % flux z element node extrap. max. value
        fm_elemextrap         double = double.empty; % flux module element node extrap. array
        fm_elemextrap_min     double = double.empty; % flux module element node extrap. min. value
        fm_elemextrap_max     double = double.empty; % flux module element node extrap. max. value
        
        % Structural analysis stresses smoothed at nodes
        sxx_nodeextrap        double = double.empty; % sigma x extrap. node smoothed stress array
        sxx_nodeextrap_min    double = double.empty; % sigma x extrap. node smoothed min. stress
        sxx_nodeextrap_max    double = double.empty; % sigma x extrap. node smoothed max. stress
        syy_nodeextrap        double = double.empty; % sigma y extrap. node smoothed stress array
        syy_nodeextrap_min    double = double.empty; % sigma y extrap. node smoothed min. stress
        syy_nodeextrap_max    double = double.empty; % sigma y extrap. node smoothed max. stress
        szz_nodeextrap        double = double.empty; % sigma z extrap. node smoothed stress array
        szz_nodeextrap_min    double = double.empty; % sigma z extrap. node smoothed min. stress
        szz_nodeextrap_max    double = double.empty; % sigma z extrap. node smoothed max. stress
        txy_nodeextrap        double = double.empty; % tau xy extrap. node smoothed stress array
        txy_nodeextrap_min    double = double.empty; % tau xy extrap. node smoothed min. stress
        txy_nodeextrap_max    double = double.empty; % tau xy extrap. node smoothed max. stress
        txz_nodeextrap        double = double.empty; % tau xz extrap. node smoothed stress array
        txz_nodeextrap_min    double = double.empty; % tau xz extrap. node smoothed min. stress
        txz_nodeextrap_max    double = double.empty; % tau xz extrap. node smoothed max. stress
        tyz_nodeextrap        double = double.empty; % tau yz extrap. node smoothed stress array
        tyz_nodeextrap_min    double = double.empty; % tau yz extrap. node smoothed min. stress
        tyz_nodeextrap_max    double = double.empty; % tau yz extrap. node smoothed max. stress
        s1_nodeextrap         double = double.empty; % sigma 1 extrap. node smoothed stress array
        s1_nodeextrap_min     double = double.empty; % sigma 1 extrap. node smoothed min. stress
        s1_nodeextrap_max     double = double.empty; % sigma 1 extrap. node smoothed max. stress
        s2_nodeextrap         double = double.empty; % sigma 2 extrap. node smoothed stress array
        s2_nodeextrap_min     double = double.empty; % sigma 2 extrap. node smoothed min. stress
        s2_nodeextrap_max     double = double.empty; % sigma 2 extrap. node smoothed max. stress
        s3_nodeextrap         double = double.empty; % sigma 3 extrap. node smoothed stress array
        s3_nodeextrap_min     double = double.empty; % sigma 3 extrap. node smoothed min. stress
        s3_nodeextrap_max     double = double.empty; % sigma 3 extrap. node smoothed max. stress
        tmax_nodeextrap       double = double.empty; % tau max. extrap. node smoothed stress array
        tmax_nodeextrap_min   double = double.empty; % tau max. extrap. node smoothed min. stress
        tmax_nodeextrap_max   double = double.empty; % tau max. extrap. node smoothed max. stress
        
        % Structural analysis internal forces smoothed at nodes
        qxz_nodeextrap        double = double.empty; % shear xz extrap. node smoothed shear array
        qxz_nodeextrap_min    double = double.empty; % shear xz extrap. node smoothed min. shear array
        qxz_nodeextrap_max    double = double.empty; % shear xz extrap. node smoothed max. shear array
        qyz_nodeextrap        double = double.empty; % shear yz extrap. node smoothed shear array
        qyz_nodeextrap_min    double = double.empty; % shear yz extrap. node smoothed min. shear array
        qyz_nodeextrap_max    double = double.empty; % shear yz extrap. node smoothed max. shear array
        mxx_nodeextrap        double = double.empty; % moment x extrap. node smoothed moments array
        mxx_nodeextrap_min    double = double.empty; % moment x extrap. node smoothed min. moments array
        mxx_nodeextrap_max    double = double.empty; % moment x extrap. node smoothed max. moments array
        myy_nodeextrap        double = double.empty; % moment y extrap. node smoothed moments array
        myy_nodeextrap_min    double = double.empty; % moment y extrap. node smoothed min. moments array
        myy_nodeextrap_max    double = double.empty; % moment y extrap. node smoothed max. moments array
        mxy_nodeextrap        double = double.empty; % moment xy extrap. node smoothed moments array
        mxy_nodeextrap_min    double = double.empty; % moment xy extrap. node smoothed min. moments array
        mxy_nodeextrap_max    double = double.empty; % moment xy extrap. node smoothed max. moments array
        m1_nodeextrap         double = double.empty; % moment 1 extrap. node smoothed moments array
        m1_nodeextrap_min     double = double.empty; % moment 1 extrap. node smoothed min. moments array
        m1_nodeextrap_max     double = double.empty; % moment 1 extrap. node smoothed max. moments array
        m2_nodeextrap         double = double.empty; % moment 2 extrap. node smoothed moments array
        m2_nodeextrap_min     double = double.empty; % moment 2 extrap. node smoothed min. moments array
        m2_nodeextrap_max     double = double.empty; % moment 2 extrap. node smoothed max. moments array
        tormax_nodeextrap     double = double.empty; % torsion max. extrap. node smoothed torsion array
        tormax_nodeextrap_min double = double.empty; % torsion max. extrap. node smoothed min. torsion array
        tormax_nodeextrap_max double = double.empty; % torsion max. extrap. node smoothed max. torsion array
        
        % Thermal analysis heat fluxes smoothed at nodes
        fxx_nodeextrap        double = double.empty; % flux x extrap. node smoothed array
        fxx_nodeextrap_min    double = double.empty; % flux x extrap. node smoothed min. value
        fxx_nodeextrap_max    double = double.empty; % flux x extrap. node smoothed max. value
        fyy_nodeextrap        double = double.empty; % flux y extrap. node smoothed array
        fyy_nodeextrap_min    double = double.empty; % flux y extrap. node smoothed min. value
        fyy_nodeextrap_max    double = double.empty; % flux y extrap. node smoothed max. value
        fzz_nodeextrap        double = double.empty; % flux z extrap. node smoothed array
        fzz_nodeextrap_min    double = double.empty; % flux z extrap. node smoothed min. value
        fzz_nodeextrap_max    double = double.empty; % flux z extrap. node smoothed max. value
        fm_nodeextrap         double = double.empty; % flux module extrap. node smoothed array
        fm_nodeextrap_min     double = double.empty; % flux module extrap. node smoothed min. value
        fm_nodeextrap_max     double = double.empty; % flux module extrap. node smoothed max. value
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
            
            if (mdl.anm.SIGMA_XX)
                this.sxx_gp         = zeros(this.maxGPts,mdl.nel);
                this.sxx_elemextrap = zeros(this.maxNen,mdl.nel);
                this.sxx_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SIGMA_YY)
                this.syy_gp         = zeros(this.maxGPts,mdl.nel);
                this.syy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.syy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SIGMA_ZZ)
                this.szz_gp         = zeros(this.maxGPts,mdl.nel);
                this.szz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.szz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.TAU_XY)
                this.txy_gp         = zeros(this.maxGPts,mdl.nel);
                this.txy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.txy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.TAU_XZ)
                this.txz_gp         = zeros(this.maxGPts,mdl.nel);
                this.txz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.txz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.TAU_YZ)
                this.tyz_gp         = zeros(this.maxGPts,mdl.nel);
                this.tyz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.tyz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SIGMA_1)
                this.s1_gp         = zeros(this.maxGPts,mdl.nel);
                this.s1x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s1y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s1_elemextrap = zeros(this.maxNen,mdl.nel);
                this.s1_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SIGMA_2)
                this.s2_gp         = zeros(this.maxGPts,mdl.nel);
                this.s2x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s2y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s2_elemextrap = zeros(this.maxNen,mdl.nel);
                this.s2_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SIGMA_3)
                this.s3_gp         = zeros(this.maxGPts,mdl.nel);
                this.s3x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s3y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.s3_elemextrap = zeros(this.maxNen,mdl.nel);
                this.s3_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.TAU_MAX)
                this.tmax_gp         = zeros(this.maxGPts,mdl.nel);
                this.tmax_elemextrap = zeros(this.maxNen,mdl.nel);
                this.tmax_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SHEAR_XZ)
                this.qxz_gp         = zeros(this.maxGPts,mdl.nel);
                this.qxz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.qxz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.SHEAR_YZ)
                this.qyz_gp         = zeros(this.maxGPts,mdl.nel);
                this.qyz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.qyz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.MOMENT_XX)
                this.mxx_gp         = zeros(this.maxGPts,mdl.nel);
                this.mxx_elemextrap = zeros(this.maxNen,mdl.nel);
                this.mxx_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.MOMENT_YY)
                this.myy_gp         = zeros(this.maxGPts,mdl.nel);
                this.myy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.myy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.MOMENT_XY)
                this.mxy_gp         = zeros(this.maxGPts,mdl.nel);
                this.mxy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.mxy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.MOMENT_1)
                this.m1_gp         = zeros(this.maxGPts,mdl.nel);
                this.m1x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.m1y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.m1_elemextrap = zeros(this.maxNen,mdl.nel);
                this.m1_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.MOMENT_2)
                this.m2_gp         = zeros(this.maxGPts,mdl.nel);
                this.m2x_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.m2y_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.m2_elemextrap = zeros(this.maxNen,mdl.nel);
                this.m2_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.TORSION_MAX)
                this.tormax_gp         = zeros(this.maxGPts,mdl.nel);
                this.tormax_elemextrap = zeros(this.maxNen,mdl.nel);
                this.tormax_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.FLUX_XX)
                this.fxx_gp         = zeros(this.maxGPts,mdl.nel);
                this.fxx_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fxx_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.FLUX_YY)
                this.fyy_gp         = zeros(this.maxGPts,mdl.nel);
                this.fyy_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fyy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.FLUX_ZZ)
                this.fzz_gp         = zeros(this.maxGPts,mdl.nel);
                this.fzz_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fzz_nodeextrap = zeros(mdl.nnp,1);
            end
            if (mdl.anm.FLUX_MOD)
                this.fm_gp         = zeros(this.maxGPts,mdl.nel);
                this.fmx_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.fmy_gp        = zeros(this.maxGPts*mdl.nel,1);
                this.fm_elemextrap = zeros(this.maxNen,mdl.nel);
                this.fm_nodeextrap = zeros(mdl.nnp,1);
            end
        end
        
        %------------------------------------------------------------------
        % Compute minimum and maximum values of obtained results.
        function setMinMaxValues(this,mdl)
            if (mdl.anm.SIGMA_XX)
                this.sxx_gp_min         = min(min(this.sxx_gp));
                this.sxx_gp_max         = max(max(this.sxx_gp));
                this.sxx_elemextrap_min = min(min(this.sxx_elemextrap));
                this.sxx_elemextrap_max = max(max(this.sxx_elemextrap));
                this.sxx_nodeextrap_min = min(min(this.sxx_nodeextrap));
                this.sxx_nodeextrap_max = max(max(this.sxx_nodeextrap));
            end
            if (mdl.anm.SIGMA_YY)
                this.syy_gp_min         = min(min(this.syy_gp));
                this.syy_gp_max         = max(max(this.syy_gp));
                this.syy_elemextrap_min = min(min(this.syy_elemextrap));
                this.syy_elemextrap_max = max(max(this.syy_elemextrap));
                this.syy_nodeextrap_min = min(min(this.syy_nodeextrap));
                this.syy_nodeextrap_max = max(max(this.syy_nodeextrap));
            end
            if (mdl.anm.SIGMA_ZZ)
                this.szz_gp_min         = min(min(this.szz_gp));
                this.szz_gp_max         = max(max(this.szz_gp));
                this.szz_elemextrap_min = min(min(this.szz_elemextrap));
                this.szz_elemextrap_max = max(max(this.szz_elemextrap));
                this.szz_nodeextrap_min = min(min(this.szz_nodeextrap));
                this.szz_nodeextrap_max = max(max(this.szz_nodeextrap));
            end
            if (mdl.anm.TAU_XY)
                this.txy_gp_min         = min(min(this.txy_gp));
                this.txy_gp_max         = max(max(this.txy_gp));
                this.txy_elemextrap_min = min(min(this.txy_elemextrap));
                this.txy_elemextrap_max = max(max(this.txy_elemextrap));
                this.txy_nodeextrap_min = min(min(this.txy_nodeextrap));
                this.txy_nodeextrap_max = max(max(this.txy_nodeextrap));
            end
            if (mdl.anm.TAU_XZ)
                this.txz_gp_min         = min(min(this.txz_gp));
                this.txz_gp_max         = max(max(this.txz_gp));
                this.txz_elemextrap_min = min(min(this.txz_elemextrap));
                this.txz_elemextrap_max = max(max(this.txz_elemextrap));
                this.txz_nodeextrap_min = min(min(this.txz_nodeextrap));
                this.txz_nodeextrap_max = max(max(this.txz_nodeextrap));
            end
            if (mdl.anm.TAU_YZ)
                this.tyz_gp_min         = min(min(this.tyz_gp));
                this.tyz_gp_max         = max(max(this.tyz_gp));
                this.tyz_elemextrap_min = min(min(this.tyz_elemextrap));
                this.tyz_elemextrap_max = max(max(this.tyz_elemextrap));
                this.tyz_nodeextrap_min = min(min(this.tyz_nodeextrap));
                this.tyz_nodeextrap_max = max(max(this.tyz_nodeextrap));
            end
            if (mdl.anm.SIGMA_1)
                this.s1_gp_min         = min(min(this.s1_gp));
                this.s1_gp_max         = max(max(this.s1_gp));
                this.s1_elemextrap_min = min(min(this.s1_elemextrap));
                this.s1_elemextrap_max = max(max(this.s1_elemextrap));
                this.s1_nodeextrap_min = min(min(this.s1_nodeextrap));
                this.s1_nodeextrap_max = max(max(this.s1_nodeextrap));
            end
            if (mdl.anm.SIGMA_2)
                this.s2_gp_min         = min(min(this.s2_gp));
                this.s2_gp_max         = max(max(this.s2_gp));
                this.s2_elemextrap_min = min(min(this.s2_elemextrap));
                this.s2_elemextrap_max = max(max(this.s2_elemextrap));
                this.s2_nodeextrap_min = min(min(this.s2_nodeextrap));
                this.s2_nodeextrap_max = max(max(this.s2_nodeextrap));
            end
            if (mdl.anm.SIGMA_3)
                this.s3_gp_min         = min(min(this.s3_gp));
                this.s3_gp_max         = max(max(this.s3_gp));
                this.s3_elemextrap_min = min(min(this.s3_elemextrap));
                this.s3_elemextrap_max = max(max(this.s3_elemextrap));
                this.s3_nodeextrap_min = min(min(this.s3_nodeextrap));
                this.s3_nodeextrap_max = max(max(this.s3_nodeextrap));
            end
            if (mdl.anm.TAU_MAX)
                this.tmax_gp_min         = min(min(this.tmax_gp));
                this.tmax_gp_max         = max(max(this.tmax_gp));
                this.tmax_elemextrap_min = min(min(this.tmax_elemextrap));
                this.tmax_elemextrap_max = max(max(this.tmax_elemextrap));
                this.tmax_nodeextrap_min = min(min(this.tmax_nodeextrap));
                this.tmax_nodeextrap_max = max(max(this.tmax_nodeextrap));
            end        
            if (mdl.anm.SHEAR_XZ)
                this.qxz_gp_min         = min(min(this.qxz_gp));
                this.qxz_gp_max         = max(max(this.qxz_gp));
                this.qxz_elemextrap_min = min(min(this.qxz_elemextrap));
                this.qxz_elemextrap_max = max(max(this.qxz_elemextrap));
                this.qxz_nodeextrap_min = min(min(this.qxz_nodeextrap));
                this.qxz_nodeextrap_max = max(max(this.qxz_nodeextrap));
            end
            if (mdl.anm.SHEAR_YZ)
                this.qyz_gp_min         = min(min(this.qyz_gp));
                this.qyz_gp_max         = max(max(this.qyz_gp));
                this.qyz_elemextrap_min = min(min(this.qyz_elemextrap));
                this.qyz_elemextrap_max = max(max(this.qyz_elemextrap));
                this.qyz_nodeextrap_min = min(min(this.qyz_nodeextrap));
                this.qyz_nodeextrap_max = max(max(this.qyz_nodeextrap));
            end
            if (mdl.anm.MOMENT_XX)
                this.mxx_gp_min         = min(min(this.mxx_gp));
                this.mxx_gp_max         = max(max(this.mxx_gp));
                this.mxx_elemextrap_min = min(min(this.mxx_elemextrap));
                this.mxx_elemextrap_max = max(max(this.mxx_elemextrap));
                this.mxx_nodeextrap_min = min(min(this.mxx_nodeextrap));
                this.mxx_nodeextrap_max = max(max(this.mxx_nodeextrap));
            end
            if (mdl.anm.MOMENT_YY)
                this.myy_gp_min         = min(min(this.myy_gp));
                this.myy_gp_max         = max(max(this.myy_gp));
                this.myy_elemextrap_min = min(min(this.myy_elemextrap));
                this.myy_elemextrap_max = max(max(this.myy_elemextrap));
                this.myy_nodeextrap_min = min(min(this.myy_nodeextrap));
                this.myy_nodeextrap_max = max(max(this.myy_nodeextrap));
            end
            if (mdl.anm.MOMENT_XY)
                this.mxy_gp_min         = min(min(this.mxy_gp));
                this.mxy_gp_max         = max(max(this.mxy_gp));
                this.mxy_elemextrap_min = min(min(this.mxy_elemextrap));
                this.mxy_elemextrap_max = max(max(this.mxy_elemextrap));
                this.mxy_nodeextrap_min = min(min(this.mxy_nodeextrap));
                this.mxy_nodeextrap_max = max(max(this.mxy_nodeextrap));
            end
            if (mdl.anm.MOMENT_1)
                this.m1_gp_min         = min(min(this.m1_gp));
                this.m1_gp_max         = max(max(this.m1_gp));
                this.m1_elemextrap_min = min(min(this.m1_elemextrap));
                this.m1_elemextrap_max = max(max(this.m1_elemextrap));
                this.m1_nodeextrap_min = min(min(this.m1_nodeextrap));
                this.m1_nodeextrap_max = max(max(this.m1_nodeextrap));
            end
            if (mdl.anm.MOMENT_2)
                this.m2_gp_min         = min(min(this.m2_gp));
                this.m2_gp_max         = max(max(this.m2_gp));
                this.m2_elemextrap_min = min(min(this.m2_elemextrap));
                this.m2_elemextrap_max = max(max(this.m2_elemextrap));
                this.m2_nodeextrap_min = min(min(this.m2_nodeextrap));
                this.m2_nodeextrap_max = max(max(this.m2_nodeextrap));
            end
            if (mdl.anm.TORSION_MAX)
                this.tormax_gp_min         = min(min(this.tormax_gp));
                this.tormax_gp_max         = max(max(this.tormax_gp));
                this.tormax_elemextrap_min = min(min(this.tormax_elemextrap));
                this.tormax_elemextrap_max = max(max(this.tormax_elemextrap));
                this.tormax_nodeextrap_min = min(min(this.tormax_nodeextrap));
                this.tormax_nodeextrap_max = max(max(this.tormax_nodeextrap));
            end
            if (mdl.anm.FLUX_XX)
                this.fxx_gp_min         = min(min(this.fxx_gp));
                this.fxx_gp_max         = max(max(this.fxx_gp));
                this.fxx_elemextrap_min = min(min(this.fxx_elemextrap));
                this.fxx_elemextrap_max = max(max(this.fxx_elemextrap));
                this.fxx_nodeextrap_min = min(min(this.fxx_nodeextrap));
                this.fxx_nodeextrap_max = max(max(this.fxx_nodeextrap));
            end
            if (mdl.anm.FLUX_YY)
                this.fyy_gp_min         = min(min(this.fyy_gp));
                this.fyy_gp_max         = max(max(this.fyy_gp));
                this.fyy_elemextrap_min = min(min(this.fyy_elemextrap));
                this.fyy_elemextrap_max = max(max(this.fyy_elemextrap));
                this.fyy_nodeextrap_min = min(min(this.fyy_nodeextrap));
                this.fyy_nodeextrap_max = max(max(this.fyy_nodeextrap));
            end
            if (mdl.anm.FLUX_ZZ)
                this.fzz_gp_min         = min(min(this.fzz_gp));
                this.fzz_gp_max         = max(max(this.fzz_gp));
                this.fzz_elemextrap_min = min(min(this.fzz_elemextrap));
                this.fzz_elemextrap_max = max(max(this.fzz_elemextrap));
                this.fzz_nodeextrap_min = min(min(this.fzz_nodeextrap));
                this.fzz_nodeextrap_max = max(max(this.fzz_nodeextrap));
            end
            if (mdl.anm.FLUX_MOD)
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
            if (mdl.anm.SIGMA_XX)
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
            if (mdl.anm.SIGMA_YY)
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
            if (mdl.anm.SIGMA_ZZ)
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
            if (mdl.anm.TAU_XY)
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
            if (mdl.anm.TAU_XZ)
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
            if (mdl.anm.TAU_YZ)
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
            if (mdl.anm.SIGMA_1)
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
            if (mdl.anm.SIGMA_2)
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
            if (mdl.anm.SIGMA_3)
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
            if (mdl.anm.TAU_MAX)
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
            %--------------------------------------------------------------
            if (mdl.anm.SHEAR_XZ)
                % Clear small values
                if (abs(this.qxz_gp_min) < this.tol && abs(this.qxz_gp_max) < this.tol)
                    this.qxz_gp_min = 0.0;
                    this.qxz_gp_max = 0.0;
                    this.qxz_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.qxz_elemextrap_min) < this.tol && abs(this.qxz_elemextrap_max) < this.tol)
                    this.qxz_elemextrap_min = 0.0;
                    this.qxz_elemextrap_max = 0.0;
                    this.qxz_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.qxz_nodeextrap_min) < this.tol && abs(this.qxz_nodeextrap_max) < this.tol)
                    this.qxz_nodeextrap_min = 0.0;
                    this.qxz_nodeextrap_max = 0.0;
                    this.qxz_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.qxz_gp_min - this.qxz_gp_max) < this.tol)
                    mean = (this.qxz_gp_min+this.qxz_gp_max)/2;
                    this.qxz_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.qxz_gp_min = this.qxz_gp_min-abs(this.qxz_gp_min/10);
                    this.qxz_gp_max = this.qxz_gp_max+abs(this.qxz_gp_max/10);
                end
                if (abs(this.qxz_elemextrap_min - this.qxz_elemextrap_max) < this.tol)
                    mean = (this.qxz_elemextrap_min+this.qxz_elemextrap_max)/2;
                    this.qxz_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.qxz_elemextrap_min = this.qxz_elemextrap_min-abs(this.qxz_elemextrap_min/10);
                    this.qxz_elemextrap_max = this.qxz_elemextrap_max+abs(this.qxz_elemextrap_max/10);
                end
                if (abs(this.qxz_nodeextrap_min - this.qxz_nodeextrap_max) < this.tol)
                    mean = (this.qxz_nodeextrap_min+this.qxz_nodeextrap_max)/2;
                    this.qxz_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.qxz_nodeextrap_min = this.qxz_nodeextrap_min-abs(this.qxz_nodeextrap_min/10);
                    this.qxz_nodeextrap_max = this.qxz_nodeextrap_max+abs(this.qxz_nodeextrap_max/10);
                end
            end
            if (mdl.anm.SHEAR_YZ)
                % Clear small values
                if (abs(this.qyz_gp_min) < this.tol && abs(this.qyz_gp_max) < this.tol)
                    this.qyz_gp_min = 0.0;
                    this.qyz_gp_max = 0.0;
                    this.qyz_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.qyz_elemextrap_min) < this.tol && abs(this.qyz_elemextrap_max) < this.tol)
                    this.qyz_elemextrap_min = 0.0;
                    this.qyz_elemextrap_max = 0.0;
                    this.qyz_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.qyz_nodeextrap_min) < this.tol && abs(this.qyz_nodeextrap_max) < this.tol)
                    this.qyz_nodeextrap_min = 0.0;
                    this.qyz_nodeextrap_max = 0.0;
                    this.qyz_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.qyz_gp_min - this.qyz_gp_max) < this.tol)
                    mean = (this.qyz_gp_min+this.qyz_gp_max)/2;
                    this.qyz_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.qyz_gp_min = this.qyz_gp_min-abs(this.qyz_gp_min/10);
                    this.qyz_gp_max = this.qyz_gp_max+abs(this.qyz_gp_max/10);
                end
                if (abs(this.qyz_elemextrap_min - this.qyz_elemextrap_max) < this.tol)
                    mean = (this.qyz_elemextrap_min+this.qyz_elemextrap_max)/2;
                    this.qyz_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.qyz_elemextrap_min = this.qyz_elemextrap_min-abs(this.qyz_elemextrap_min/10);
                    this.qyz_elemextrap_max = this.qyz_elemextrap_max+abs(this.qyz_elemextrap_max/10);
                end
                if (abs(this.qyz_nodeextrap_min - this.qyz_nodeextrap_max) < this.tol)
                    mean = (this.qyz_nodeextrap_min+this.qyz_nodeextrap_max)/2;
                    this.qyz_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.qyz_nodeextrap_min = this.qyz_nodeextrap_min-abs(this.qyz_nodeextrap_min/10);
                    this.qyz_nodeextrap_max = this.qyz_nodeextrap_max+abs(this.qyz_nodeextrap_max/10);
                end
            end
            if (mdl.anm.MOMENT_XX)
                % Clear small values
                if (abs(this.mxx_gp_min) < this.tol && abs(this.mxx_gp_max) < this.tol)
                    this.mxx_gp_min = 0.0;
                    this.mxx_gp_max = 0.0;
                    this.mxx_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.mxx_elemextrap_min) < this.tol && abs(this.mxx_elemextrap_max) < this.tol)
                    this.mxx_elemextrap_min = 0.0;
                    this.mxx_elemextrap_max = 0.0;
                    this.mxx_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.mxx_nodeextrap_min) < this.tol && abs(this.mxx_nodeextrap_max) < this.tol)
                    this.mxx_nodeextrap_min = 0.0;
                    this.mxx_nodeextrap_max = 0.0;
                    this.mxx_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.mxx_gp_min - this.mxx_gp_max) < this.tol)
                    mean = (this.mxx_gp_min+this.mxx_gp_max)/2;
                    this.mxx_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.mxx_gp_min = this.mxx_gp_min-abs(this.mxx_gp_min/10);
                    this.mxx_gp_max = this.mxx_gp_max+abs(this.mxx_gp_max/10);
                end
                if (abs(this.mxx_elemextrap_min - this.mxx_elemextrap_max) < this.tol)
                    mean = (this.mxx_elemextrap_min+this.mxx_elemextrap_max)/2;
                    this.mxx_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.mxx_elemextrap_min = this.mxx_elemextrap_min-abs(this.mxx_elemextrap_min/10);
                    this.mxx_elemextrap_max = this.mxx_elemextrap_max+abs(this.mxx_elemextrap_max/10);
                end
                if (abs(this.mxx_nodeextrap_min - this.mxx_nodeextrap_max) < this.tol)
                    mean = (this.mxx_nodeextrap_min+this.mxx_nodeextrap_max)/2;
                    this.mxx_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.mxx_nodeextrap_min = this.mxx_nodeextrap_min-abs(this.mxx_nodeextrap_min/10);
                    this.mxx_nodeextrap_max = this.mxx_nodeextrap_max+abs(this.mxx_nodeextrap_max/10);
                end
            end
            if (mdl.anm.MOMENT_YY)
                % Clear small values
                if (abs(this.myy_gp_min) < this.tol && abs(this.myy_gp_max) < this.tol)
                    this.myy_gp_min = 0.0;
                    this.myy_gp_max = 0.0;
                    this.myy_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.myy_elemextrap_min) < this.tol && abs(this.myy_elemextrap_max) < this.tol)
                    this.myy_elemextrap_min = 0.0;
                    this.myy_elemextrap_max = 0.0;
                    this.myy_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.myy_nodeextrap_min) < this.tol && abs(this.myy_nodeextrap_max) < this.tol)
                    this.myy_nodeextrap_min = 0.0;
                    this.myy_nodeextrap_max = 0.0;
                    this.myy_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.myy_gp_min - this.myy_gp_max) < this.tol)
                    mean = (this.myy_gp_min+this.myy_gp_max)/2;
                    this.myy_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.myy_gp_min = this.myy_gp_min-abs(this.myy_gp_min/10);
                    this.myy_gp_max = this.myy_gp_max+abs(this.myy_gp_max/10);
                end
                if (abs(this.myy_elemextrap_min - this.myy_elemextrap_max) < this.tol)
                    mean = (this.myy_elemextrap_min+this.myy_elemextrap_max)/2;
                    this.myy_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.myy_elemextrap_min = this.myy_elemextrap_min-abs(this.myy_elemextrap_min/10);
                    this.myy_elemextrap_max = this.myy_elemextrap_max+abs(this.myy_elemextrap_max/10);
                end
                if (abs(this.myy_nodeextrap_min - this.myy_nodeextrap_max) < this.tol)
                    mean = (this.myy_nodeextrap_min+this.myy_nodeextrap_max)/2;
                    this.myy_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.myy_nodeextrap_min = this.myy_nodeextrap_min-abs(this.myy_nodeextrap_min/10);
                    this.myy_nodeextrap_max = this.myy_nodeextrap_max+abs(this.myy_nodeextrap_max/10);
                end
            end
            if (mdl.anm.MOMENT_XY)
                % Clear small values
                if (abs(this.mxy_gp_min) < this.tol && abs(this.mxy_gp_max) < this.tol)
                    this.mxy_gp_min = 0.0;
                    this.mxy_gp_max = 0.0;
                    this.mxy_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.mxy_elemextrap_min) < this.tol && abs(this.mxy_elemextrap_max) < this.tol)
                    this.mxy_elemextrap_min = 0.0;
                    this.mxy_elemextrap_max = 0.0;
                    this.mxy_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.mxy_nodeextrap_min) < this.tol && abs(this.mxy_nodeextrap_max) < this.tol)
                    this.mxy_nodeextrap_min = 0.0;
                    this.mxy_nodeextrap_max = 0.0;
                    this.mxy_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.mxy_gp_min - this.mxy_gp_max) < this.tol)
                    mean = (this.mxy_gp_min+this.mxy_gp_max)/2;
                    this.mxy_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.mxy_gp_min = this.mxy_gp_min-abs(this.mxy_gp_min/10);
                    this.mxy_gp_max = this.mxy_gp_max+abs(this.mxy_gp_max/10);
                end
                if (abs(this.mxy_elemextrap_min - this.mxy_elemextrap_max) < this.tol)
                    mean = (this.mxy_elemextrap_min+this.mxy_elemextrap_max)/2;
                    this.mxy_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.mxy_elemextrap_min = this.mxy_elemextrap_min-abs(this.mxy_elemextrap_min/10);
                    this.mxy_elemextrap_max = this.mxy_elemextrap_max+abs(this.mxy_elemextrap_max/10);
                end
                if (abs(this.mxy_nodeextrap_min - this.mxy_nodeextrap_max) < this.tol)
                    mean = (this.mxy_nodeextrap_min+this.mxy_nodeextrap_max)/2;
                    this.mxy_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.mxy_nodeextrap_min = this.mxy_nodeextrap_min-abs(this.mxy_nodeextrap_min/10);
                    this.mxy_nodeextrap_max = this.mxy_nodeextrap_max+abs(this.mxy_nodeextrap_max/10);
                end
            end
            if (mdl.anm.MOMENT_1)
                % Clear small values
                if (abs(this.m1_gp_min) < this.tol && abs(this.m1_gp_max) < this.tol)
                    this.m1_gp_min = 0.0;
                    this.m1_gp_max = 0.0;
                    this.m1_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.m1_elemextrap_min) < this.tol && abs(this.m1_elemextrap_max) < this.tol)
                    this.m1_elemextrap_min = 0.0;
                    this.m1_elemextrap_max = 0.0;
                    this.m1_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.m1_nodeextrap_min) < this.tol && abs(this.m1_nodeextrap_max) < this.tol)
                    this.m1_nodeextrap_min = 0.0;
                    this.m1_nodeextrap_max = 0.0;
                    this.m1_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.m1_gp_min - this.m1_gp_max) < this.tol)
                    mean = (this.m1_gp_min+this.m1_gp_max)/2;
                    this.m1_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.m1_gp_min = this.m1_gp_min-abs(this.m1_gp_min/10);
                    this.m1_gp_max = this.m1_gp_max+abs(this.m1_gp_max/10);
                end
                if (abs(this.m1_elemextrap_min - this.m1_elemextrap_max) < this.tol)
                    mean = (this.m1_elemextrap_min+this.m1_elemextrap_max)/2;
                    this.m1_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.m1_elemextrap_min = this.m1_elemextrap_min-abs(this.m1_elemextrap_min/10);
                    this.m1_elemextrap_max = this.m1_elemextrap_max+abs(this.m1_elemextrap_max/10);
                end
                if (abs(this.m1_nodeextrap_min - this.m1_nodeextrap_max) < this.tol)
                    mean = (this.m1_nodeextrap_min+this.m1_nodeextrap_max)/2;
                    this.m1_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.m1_nodeextrap_min = this.m1_nodeextrap_min-abs(this.m1_nodeextrap_min/10);
                    this.m1_nodeextrap_max = this.m1_nodeextrap_max+abs(this.m1_nodeextrap_max/10);
                end
            end
            if (mdl.anm.MOMENT_2)
                % Clear small values
                if (abs(this.m2_gp_min) < this.tol && abs(this.m2_gp_max) < this.tol)
                    this.m2_gp_min = 0.0;
                    this.m2_gp_max = 0.0;
                    this.m2_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.m2_elemextrap_min) < this.tol && abs(this.m2_elemextrap_max) < this.tol)
                    this.m2_elemextrap_min = 0.0;
                    this.m2_elemextrap_max = 0.0;
                    this.m2_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.m2_nodeextrap_min) < this.tol && abs(this.m2_nodeextrap_max) < this.tol)
                    this.m2_nodeextrap_min = 0.0;
                    this.m2_nodeextrap_max = 0.0;
                    this.m2_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.m2_gp_min - this.m2_gp_max) < this.tol)
                    mean = (this.m2_gp_min+this.m2_gp_max)/2;
                    this.m2_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.m2_gp_min = this.m2_gp_min-abs(this.m2_gp_min/10);
                    this.m2_gp_max = this.m2_gp_max+abs(this.m2_gp_max/10);
                end
                if (abs(this.m2_elemextrap_min - this.m2_elemextrap_max) < this.tol)
                    mean = (this.m2_elemextrap_min+this.m2_elemextrap_max)/2;
                    this.m2_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.m2_elemextrap_min = this.m2_elemextrap_min-abs(this.m2_elemextrap_min/10);
                    this.m2_elemextrap_max = this.m2_elemextrap_max+abs(this.m2_elemextrap_max/10);
                end
                if (abs(this.m2_nodeextrap_min - this.m2_nodeextrap_max) < this.tol)
                    mean = (this.m2_nodeextrap_min+this.m2_nodeextrap_max)/2;
                    this.m2_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.m2_nodeextrap_min = this.m2_nodeextrap_min-abs(this.m2_nodeextrap_min/10);
                    this.m2_nodeextrap_max = this.m2_nodeextrap_max+abs(this.m2_nodeextrap_max/10);
                end
                
            end
            if (mdl.anm.TORSION_MAX)
                % Clear small values
                if (abs(this.tormax_gp_min) < this.tol && abs(this.tormax_gp_max) < this.tol)
                    this.tormax_gp_min = 0.0;
                    this.tormax_gp_max = 0.0;
                    this.tormax_gp     = zeros(this.maxGPts,mdl.nel);
                end
                if (abs(this.tormax_elemextrap_min) < this.tol && abs(this.tormax_elemextrap_max) < this.tol)
                    this.tormax_elemextrap_min = 0.0;
                    this.tormax_elemextrap_max = 0.0;
                    this.tormax_elemextrap     = zeros(this.maxNen,mdl.nel);
                end
                if (abs(this.tormax_nodeextrap_min) < this.tol && abs(this.tormax_nodeextrap_max) < this.tol)
                    this.tormax_nodeextrap_min = 0.0;
                    this.tormax_nodeextrap_max = 0.0;
                    this.tormax_nodeextrap     = zeros(mdl.nnp,1);
                end
                % Smooth constant values
                if (abs(this.tormax_gp_min - this.tormax_gp_max) < this.tol)
                    mean = (this.tormax_gp_min+this.tormax_gp_max)/2;
                    this.tormax_gp     = mean*ones(this.maxGPts,mdl.nel);
                    this.tormax_gp_min = this.tormax_gp_min-abs(this.tormax_gp_min/10);
                    this.tormax_gp_max = this.tormax_gp_max+abs(this.tormax_gp_max/10);
                end
                if (abs(this.tormax_elemextrap_min - this.tormax_elemextrap_max) < this.tol)
                    mean = (this.tormax_elemextrap_min+this.tormax_elemextrap_max)/2;
                    this.tormax_elemextrap     = mean*ones(this.maxNen,mdl.nel);
                    this.tormax_elemextrap_min = this.tormax_elemextrap_min-abs(this.tormax_elemextrap_min/10);
                    this.tormax_elemextrap_max = this.tormax_elemextrap_max+abs(this.tormax_elemextrap_max/10);
                end
                if (abs(this.tormax_nodeextrap_min - this.tormax_nodeextrap_max) < this.tol)
                    mean = (this.tormax_nodeextrap_min+this.tormax_nodeextrap_max)/2;
                    this.tormax_nodeextrap     = mean*ones(mdl.nnp,1);
                    this.tormax_nodeextrap_min = this.tormax_nodeextrap_min-abs(this.tormax_nodeextrap_min/10);
                    this.tormax_nodeextrap_max = this.tormax_nodeextrap_max+abs(this.tormax_nodeextrap_max/10);
                end
                
            end
            %--------------------------------------------------------------
            if (mdl.anm.FLUX_XX)
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
            if (mdl.anm.FLUX_YY)
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
            if (mdl.anm.FLUX_ZZ)
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
            if (mdl.anm.FLUX_MOD)
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
