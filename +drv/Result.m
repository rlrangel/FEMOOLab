%% Result Class
%
%% Description
%
% This class defines a result object in the FEMOOLab program.
% A Result object is responsible for storing the analysis results and
% provide them in text or graphical formats.
%
classdef Result < handle
    %% Flags for response display options
    properties (SetAccess = public, GetAccess = public)
        eid    = logical(false);  % Plot element numbers
        nid    = logical(false);  % Plot node numbers
        scl    = double(1.0);     % Scale factor for deformed mesh
        dx     = logical(false);  % Plot contour of displacements in X direction
        dy     = logical(false);  % Plot contour of displacements in Y direction
        dz     = logical(false);  % Plot contour of displacements in Z directiont
        rx     = logical(false);  % Plot contour of rotations about X axis
        ry     = logical(false);  % Plot contour of rotations about Y axis
        rz     = logical(false);  % Plot contour of rotations about Z axis
        temp   = logical(false);  % Plot contour of temperature
        smooth = logical(false);  % Smooth element results at common nodes
        tol    = double(1e-5);    % Tolerance for cleaning small result values and differences
        sxx    = logical(false);  % Plot contour of normal stresses in X direction
        syy    = logical(false);  % Plot contour of normal stresses in Y direction
        szz    = logical(false);  % Plot contour of normal stresses in Z direction
        txy    = logical(false);  % Plot contour of XY shear stresses
        txz    = logical(false);  % Plot contour of XZ shear stresses
        tyz    = logical(false);  % Plot contour of YZ shear stresses
        s1     = logical(false);  % Plot contour of principal stresses 1
        s2     = logical(false);  % Plot contour of principal stresses 2
        s3     = logical(false);  % Plot contour of principal stresses 3
        taumax = logical(false);  % Plot contour of maximum shear stresses
        mxx    = logical(false);  % Plot contour of moment about X direction
        myy    = logical(false);  % Plot contour of moment about Y direction
        mxy    = logical(false);  % Plot contour of torsion moment
        qxz    = logical(false);  % Plot contour of XZ shear force
        qyz    = logical(false);  % Plot contour of YZ shear force
        m1     = logical(false);  % Plot contour of principal moment 1
        m2     = logical(false);  % Plot contour of principal moment 2
        tormax = logical(false);  % Plot contour of maximum torsion
        fxx    = logical(false);  % Plot contour of heat fluxes in X direction
        fyy    = logical(false);  % Plot contour of heat fluxes in Y direction
        fzz    = logical(false);  % Plot contour of heat fluxes in Z direction
        fp     = logical(false);  % Plot contour of fluxes directions
    end
    
    %% Constant values for contour types
    properties (Constant = true, Access = public)
        SXX_ELEMEXTRAP  = int32(1);
        SYY_ELEMEXTRAP  = int32(2);
        TXY_ELEMEXTRAP  = int32(3);
        S1_ELEMEXTRAP   = int32(4);
        S2_ELEMEXTRAP   = int32(5);
        TMAX_ELEMEXTRAP = int32(6);
        SXX_NODEEXTRAP  = int32(7);
        SYY_NODEEXTRAP  = int32(8);
        TXY_NODEEXTRAP  = int32(9);
        S1_NODEEXTRAP   = int32(10);
        S2_NODEEXTRAP   = int32(11);
        TMAX_NODEEXTRAP = int32(12);
        TEMP_NODE       = int32(13);
        FXX_ELEMEXTRAP  = int32(14);
        FYY_ELEMEXTRAP  = int32(15);
        FP_ELEMEXTRAP   = int32(16);
        FXX_NODEEXTRAP  = int32(17);
        FYY_NODEEXTRAP  = int32(18);
        FP_NODEEXTRAP   = int32(19);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Results from equilibrium system
        U double = double.empty;                   % global displacement vector
        
        % Gauss point results
        ngp                 int32  = int32.empty;  % vector of number of element gauss pts
        x_gp                double = double.empty; % vector of gauss point x coordinates
        y_gp                double = double.empty; % vector of gauss point y coordinates
        sxx_gp              double = double.empty; % sigma x gauss points stress array
        sxx_gp_min          double = double.empty; % sigma x gauss points min. stress
        sxx_gp_max          double = double.empty; % sigma x gauss points max. stress
        syy_gp              double = double.empty; % sigma y gauss points stress array
        syy_gp_min          double = double.empty; % sigma y gauss points min. stress
        syy_gp_max          double = double.empty; % sigma y gauss points max. stress
        txy_gp              double = double.empty; % tau xy gauss points stress array
        txy_gp_min          double = double.empty; % tau xy gauss points min. stress
        txy_gp_max          double = double.empty; % tau xy gauss points max. stress
        s1_gp               double = double.empty; % sigma 1 gauss points stress array
        s1_gp_min           double = double.empty; % sigma 1 gauss points min. stress
        s1_gp_max           double = double.empty; % sigma 1 gauss points max. stress
        s2_gp               double = double.empty; % sigma 2 gauss points stress array
        s2_gp_min           double = double.empty; % sigma 2 gauss points min. stress
        s2_gp_max           double = double.empty; % sigma 2 gauss points max. stress
        tmax_gp             double = double.empty; % tau max. gauss points stress array
        tmax_gp_min         double = double.empty; % tau max. gauss points min. stress
        tmax_gp_max         double = double.empty; % tau max. gauss points max. stress
        s1x_gp              double = double.empty; % vector of s1 vector gauss x components
        s1y_gp              double = double.empty; % vector of s1 vector gauss y components
        s2x_gp              double = double.empty; % vector of s2 vector gauss x components
        s2y_gp              double = double.empty; % vector of s2 vector gauss y components
        fxx_gp              double = double.empty;
        fxx_gp_min          double = double.empty;
        fxx_gp_max          double = double.empty;
        fyy_gp              double = double.empty;
        fyy_gp_min          double = double.empty;
        fyy_gp_max          double = double.empty;
        fp_gp               double = double.empty;
        fp_gp_min           double = double.empty;
        fp_gp_max           double = double.empty;
        fpx_gp              double = double.empty;
        fpy_gp              double = double.empty;
        
        % Element node results
        sxx_elemextrap      double = double.empty; % sigma x element node extrap. stress array
        sxx_elemextrap_min  double = double.empty; % sigma x element node extrap. min. stress
        sxx_elemextrap_max  double = double.empty; % sigma x element node extrap. max. stress
        syy_elemextrap      double = double.empty; % sigma y element node extrap. stress array
        syy_elemextrap_min  double = double.empty; % sigma y element node extrap. min. stress
        syy_elemextrap_max  double = double.empty; % sigma y element node extrap. max. stress
        txy_elemextrap      double = double.empty; % tau xy element node extrap. stress array
        txy_elemextrap_min  double = double.empty; % tau xy element node extrap. min. stress
        txy_elemextrap_max  double = double.empty; % tau xy element node extrap. max. stress
        s1_elemextrap       double = double.empty; % sigma 1 element node extrap. stress array
        s1_elemextrap_min   double = double.empty; % sigma 1 element node extrap. min. stress
        s1_elemextrap_max   double = double.empty; % sigma 1 element node extrap. max. stress
        s2_elemextrap       double = double.empty; % sigma 2 element node extrap. stress array
        s2_elemextrap_min   double = double.empty; % sigma 2 element node extrap. min. stress
        s2_elemextrap_max   double = double.empty; % sigma 2 element node extrap. max. stress
        tmax_elemextrap     double = double.empty; % tau max. element node extrap. stress array
        tmax_elemextrap_min double = double.empty; % tau max. element node extrap. min stress
        tmax_elemextrap_max double = double.empty; % tau max. element node extrap. max stress
        fxx_elemextrap      double = double.empty;
        fxx_elemextrap_min  double = double.empty;
        fxx_elemextrap_max  double = double.empty;
        fyy_elemextrap      double = double.empty;
        fyy_elemextrap_min  double = double.empty;
        fyy_elemextrap_max  double = double.empty;
        fp_elemextrap       double = double.empty;
        fp_elemextrap_min   double = double.empty;
        fp_elemextrap_max   double = double.empty;
        
        % Global node results
        sxx_nodeextrap      double = double.empty; % sigma x extrap. node smoothed stress array
        sxx_nodeextrap_min  double = double.empty; % sigma x extrap. node smoothed min. stress
        sxx_nodeextrap_max  double = double.empty; % sigma x extrap. node smoothed max. stress
        syy_nodeextrap      double = double.empty; % sigma y extrap. node smoothed stress array
        syy_nodeextrap_min  double = double.empty; % sigma y extrap. node smoothed min. stress
        syy_nodeextrap_max  double = double.empty; % sigma y extrap. node smoothed max. stress
        txy_nodeextrap      double = double.empty; % tau xy extrap. node smoothed stress array
        txy_nodeextrap_min  double = double.empty; % tau xy extrap. node smoothed min. stress
        txy_nodeextrap_max  double = double.empty; % tau xy extrap. node smoothed max. stress
        s1_nodeextrap       double = double.empty; % sigma 1 extrap. node smoothed stress array
        s1_nodeextrap_min   double = double.empty; % sigma 1 extrap. node smoothed min. stress
        s1_nodeextrap_max   double = double.empty; % sigma 1 extrap. node smoothed max. stress
        s2_nodeextrap       double = double.empty; % sigma 2 extrap. node smoothed stress array
        s2_nodeextrap_min   double = double.empty; % sigma 2 extrap. node smoothed min. stress
        s2_nodeextrap_max   double = double.empty; % sigma 2 extrap. node smoothed max. stress
        tmax_nodeextrap     double = double.empty; % tau max. extrap. node smoothed stress array
        tmax_nodeextrap_min double = double.empty; % tau max. extrap. node smoothed min. stress
        tmax_nodeextrap_max double = double.empty; % tau max. extrap. node smoothed max. stress
        fxx_nodeextrap      double = double.empty;
        fxx_nodeextrap_min  double = double.empty;
        fxx_nodeextrap_max  double = double.empty;
        fyy_nodeextrap      double = double.empty;
        fyy_nodeextrap_min  double = double.empty;
        fyy_nodeextrap_max  double = double.empty;
        fp_nodeextrap       double = double.empty;
        fp_nodeextrap_min   double = double.empty;
        fp_nodeextrap_max   double = double.empty;
        
        % Handles to plot figures
        fig_deform = [];    % figure for mesh and deformed mesh plot
        fig_strbar = [];    % figure for stress bar response plots
        fig_sxx    = [];    % figure for sigma x plot
        fig_syy    = [];    % figure for sigma y plot
        fig_txy    = [];    % figure for tau xy plot
        fig_s1     = []     % figure for sigma 1 plot
        fig_s2     = []     % figure for sigma 2 plot
        fig_tmax   = []     % figure for tau max. plot
        fig_temp   = []     % figure for temperature field
        fig_fxx    = [];    % 
        fig_fyy    = [];    % 
        fig_fp     = [];    % 
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Result()
            return;
        end
    end
    
    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Clear numerical garbage from analysis. (TEMPORARY SIMPLIFICATION !!)
        function clearSmallValuesInplane(this,mdl)
            if (abs(this.sxx_gp_min) < this.tol && abs(this.sxx_gp_max) < this.tol)
                this.sxx_gp_min = 0.0;
                this.sxx_gp_max = 0.0;
                this.sxx_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.syy_gp_min) < this.tol && abs(this.syy_gp_max) < this.tol)
                this.syy_gp_min = 0.0;
                this.syy_gp_max = 0.0;
                this.syy_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.txy_gp_min) < this.tol && abs(this.txy_gp_max) < this.tol)
                this.txy_gp_min = 0.0;
                this.txy_gp_max = 0.0;
                this.txy_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.s1_gp_min) < this.tol && abs(this.s1_gp_max) < this.tol)
                this.s1_gp_min = 0.0;
                this.s1_gp_max = 0.0;
                this.s1_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.s2_gp_min) < this.tol && abs(this.s2_gp_max) < this.tol)
                this.s2_gp_min = 0.0;
                this.s2_gp_max = 0.0;
                this.s2_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.tmax_gp_min) < this.tol && abs(this.tmax_gp_max) < this.tol)
                this.tmax_gp_min = 0.0;
                this.tmax_gp_max = 0.0;
                this.tmax_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.sxx_elemextrap_min) < this.tol && abs(this.sxx_elemextrap_max) < this.tol)
                this.sxx_elemextrap_min = 0.0;
                this.sxx_elemextrap_max = 0.0;
                this.sxx_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.syy_elemextrap_min) < this.tol && abs(this.syy_elemextrap_max) < this.tol)
                this.syy_elemextrap_min = 0.0;
                this.syy_elemextrap_max = 0.0;
                this.syy_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.txy_elemextrap_min) < this.tol && abs(this.txy_elemextrap_max) < this.tol)
                this.txy_elemextrap_min = 0.0;
                this.txy_elemextrap_max = 0.0;
                this.txy_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.s1_elemextrap_min) < this.tol && abs(this.s1_elemextrap_max) < this.tol)
                this.s1_elemextrap_min = 0.0;
                this.s1_elemextrap_max = 0.0;
                this.s1_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.s2_elemextrap_min) < this.tol && abs(this.s2_elemextrap_max) < this.tol)
                this.s2_elemextrap_min = 0.0;
                this.s2_elemextrap_max = 0.0;
                this.s2_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.tmax_elemextrap_min) < this.tol && abs(this.tmax_elemextrap_max) < this.tol)
                this.tmax_elemextrap_min = 0.0;
                this.tmax_elemextrap_max = 0.0;
                this.tmax_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
        end
        
        function clearSmallValuesInplaneConduction(this,mdl)
            if (abs(this.fxx_gp_min) < this.tol && abs(this.fxx_gp_max) < this.tol)
                this.fxx_gp_min = 0.0;
                this.fxx_gp_max = 0.0;
                this.fxx_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.fyy_gp_min) < this.tol && abs(this.fyy_gp_max) < this.tol)
                this.fyy_gp_min = 0.0;
                this.fyy_gp_max = 0.0;
                this.fyy_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.fp_gp_min) < this.tol && abs(this.fp_gp_max) < this.tol)
                this.fp_gp_min = 0.0;
                this.fp_gp_max = 0.0;
                this.fp_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(this.fxx_elemextrap_min) < this.tol && abs(this.fxx_elemextrap_max) < this.tol)
                this.fxx_elemextrap_min = 0.0;
                this.fxx_elemextrap_max = 0.0;
                this.fxx_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.fyy_elemextrap_min) < this.tol && abs(this.fyy_elemextrap_max) < this.tol)
                this.fyy_elemextrap_min = 0.0;
                this.fyy_elemextrap_max = 0.0;
                this.fyy_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.fp_elemextrap_min) < this.tol && abs(this.fp_elemextrap_max) < this.tol)
                this.fp_elemextrap_min = 0.0;
                this.fp_elemextrap_max = 0.0;
                this.fp_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(this.fxx_nodeextrap_min) < this.tol && abs(this.fxx_nodeextrap_max) < this.tol)
                this.fxx_nodeextrap_min = 0.0;
                this.fxx_nodeextrap_max = 0.0;
                this.fxx_nodeextrap = zeros(mdl.nnp,1);
            end
            if (abs(this.fyy_nodeextrap_min) < this.tol && abs(this.fyy_nodeextrap_max) < this.tol)
                this.fyy_nodeextrap_min = 0.0;
                this.fyy_nodeextrap_max = 0.0;
                this.fyy_nodeextrap = zeros(mdl.nnp,1);
            end
            if (abs(this.fp_nodeextrap_min) < this.tol && abs(this.fp_nodeextrap_max) < this.tol)
                this.fp_nodeextrap_min = 0.0;
                this.fp_nodeextrap_max = 0.0;
                this.fp_nodeextrap = zeros(mdl.nnp,1);
            end
            
            % Smooth constant results
            if (abs(this.fxx_gp_min - this.fxx_gp_max) < this.tol)
                mean = (this.fxx_gp_min+this.fxx_gp_max)/2;
                this.fxx_gp = mean*ones(mdl.elems(1).gstress_npts,mdl.nel);
                this.fxx_gp_min = this.fxx_gp_min-abs(this.fxx_gp_min/10);
                this.fxx_gp_max = this.fxx_gp_max+abs(this.fxx_gp_max/10);
            end
            if (abs(this.fyy_gp_min - this.fyy_gp_max) < this.tol)
                mean = (this.fyy_gp_min+this.fyy_gp_max)/2;
                this.fxx_gp = mean*ones(mdl.elems(1).gstress_npts,mdl.nel);
                this.fyy_gp_min = this.fyy_gp_min-abs(this.fyy_gp_min/10);
                this.fyy_gp_max = this.fyy_gp_max+abs(this.fyy_gp_max/10);
            end
            if (abs(this.fp_gp_min - this.fp_gp_max) < this.tol)
                mean = (this.fp_gp_min+this.fp_gp_max)/2;
                this.fxx_gp = mean*ones(mdl.elems(1).gstress_npts,mdl.nel);
                this.fp_gp_min = this.fp_gp_min-abs(this.fp_gp_min/10);
                this.fp_gp_max = this.fp_gp_max+abs(this.fp_gp_max/10);
            end
            if (abs(this.fxx_elemextrap_min - this.fxx_elemextrap_max) < this.tol)
                mean = (this.fxx_elemextrap_min+this.fxx_elemextrap_max)/2;
                this.fxx_gp = mean*ones(mdl.elems(1).shape.nen,mdl.nel);
                this.fxx_elemextrap_min = this.fxx_elemextrap_min-abs(this.fxx_elemextrap_min/10);
                this.fxx_elemextrap_max = this.fxx_elemextrap_max+abs(this.fxx_elemextrap_max/10);
            end
            if (abs(this.fyy_elemextrap_min - this.fyy_elemextrap_max) < this.tol)
                mean = (this.fyy_elemextrap_min+this.fyy_elemextrap_max)/2;
                this.fxx_gp = mean*ones(mdl.elems(1).shape.nen,mdl.nel);
                this.fyy_elemextrap_min = this.fyy_elemextrap_min-abs(this.fyy_elemextrap_min/10);
                this.fyy_elemextrap_max = this.fyy_elemextrap_max+abs(this.fyy_elemextrap_max/10);
            end
            if (abs(this.fp_elemextrap_min - this.fp_elemextrap_max) < this.tol)
                mean = (this.fp_elemextrap_min+this.fp_elemextrap_max)/2;
                this.fxx_gp = mean*ones(mdl.elems(1).shape.nen,mdl.nel);
                this.fp_elemextrap_min = this.fp_elemextrap_min-abs(this.fp_elemextrap_min/10);
                this.fp_elemextrap_max = this.fp_elemextrap_max+abs(this.fp_elemextrap_max/10);
            end
            if (abs(this.fxx_nodeextrap_min - this.fxx_nodeextrap_max) < this.tol)
                mean = (this.fxx_nodeextrap_min+this.fxx_nodeextrap_max)/2;
                this.fxx_gp = mean*ones(mdl.nnp,1);
                this.fxx_nodeextrap_min = this.fxx_nodeextrap_min-abs(this.fxx_nodeextrap_min/10);
                this.fxx_nodeextrap_max = this.fxx_nodeextrap_max+abs(this.fxx_nodeextrap_max/10);
            end
            if (abs(this.fyy_nodeextrap_min - this.fyy_nodeextrap_max) < this.tol)
                mean = (this.fyy_nodeextrap_min+this.fyy_nodeextrap_max)/2;
                this.fxx_gp = mean*ones(mdl.nnp,1);
                this.fyy_nodeextrap_min = this.fyy_nodeextrap_min-abs(this.fyy_nodeextrap_min/10);
                this.fyy_nodeextrap_max = this.fyy_nodeextrap_max+abs(this.fyy_nodeextrap_max/10);
            end
            if (abs(this.fp_nodeextrap_min - this.fp_nodeextrap_max) < this.tol)
                mean = (this.fp_nodeextrap_min+this.fp_nodeextrap_max)/2;
                this.fxx_gp = mean*ones(mdl.nnp,1);
                this.fp_nodeextrap_min = this.fp_nodeextrap_min-abs(this.fp_nodeextrap_min/10);
                this.fp_nodeextrap_max = this.fp_nodeextrap_max+abs(this.fp_nodeextrap_max/10);
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot analysis results.
        function plotInplane(this,mdl)
            % Assemble vectors of nodal coordinates
            x = zeros(mdl.nnp,1);
            y = zeros(mdl.nnp,1);
            for i = 1:mdl.nnp
                x(i) = mdl.nodes(i).coord(1);
                y(i) = mdl.nodes(i).coord(2);
            end
            
            % Setup bounding box for displaying results
            min_x = min(x);
            max_x = max(x);
            min_y = min(y);
            max_y = max(y);
            size_x = max_x-min_x;
            size_y = max_y-min_y;
            cx = min_x + size_x*0.5;
            cy = min_y + size_y*0.5;
            plot_xmin = cx - size_x * 0.55;
            plot_xmax = cx + size_x * 0.55;
            plot_ymin = cy - size_y * 0.55;
            plot_ymax = cy + size_y * 0.55;
            
            % Create figures (windows) for displaying results
            deform_fac = this.createFigs(mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax);
            
            % Display deformed mesh
            figure(this.fig_deform);
            this.plotMesh(mdl,x,y,'k');
            this.plotDeformMesh(mdl,x,y,deform_fac,'b');
            
            % Display principal stress vectors.
            figure(this.fig_strbar);
            this.plotMesh(mdl,x,y,'k');
            quiver(this.x_gp,this.y_gp,this.s1x_gp,this.s1y_gp,'r');
            quiver(this.x_gp,this.y_gp,this.s2x_gp,this.s2y_gp,'b');
            
            % Display stress results
            if this.sxx
                figure(this.fig_sxx);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.SXX_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.SXX_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.syy
                figure(this.fig_syy);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.SYY_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.SYY_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.txy
                figure(this.fig_txy);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.TXY_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.TXY_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.s1
                figure(this.fig_s1);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.S1_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.S1_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.s2
                figure(this.fig_s2);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.S2_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.S2_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.taumax
                figure(this.fig_tmax);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.TMAX_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.TMAX_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot heat conduction analysis results.
        function plotInplaneConduction(this,mdl)
            % Assemble vectors of nodal coordinates
            x = zeros(mdl.nnp,1);
            y = zeros(mdl.nnp,1);
            for i = 1:mdl.nnp
                x(i) = mdl.nodes(i).coord(1);
                y(i) = mdl.nodes(i).coord(2);
            end
            
            % Setup bounding box for displaying results
            min_x = min(x);
            max_x = max(x);
            min_y = min(y);
            max_y = max(y);
            size_x = max_x-min_x;
            size_y = max_y-min_y;
            cx = min_x + size_x*0.5;
            cy = min_y + size_y*0.5;
            plot_xmin = cx - size_x * 0.55;
            plot_xmax = cx + size_x * 0.55;
            plot_ymin = cy - size_y * 0.55;
            plot_ymax = cy + size_y * 0.55;
            
            % Create figures (windows) for displaying results
            this.createFigs(mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax);
            
            % Display temperature results
            if this.temp
                figure(this.fig_temp);
                this.plotNodeContourInplane(mdl,x,y,drv.Result.TEMP_NODE);
                this.plotMesh(mdl,x,y,'k');
            end
            
            % Display principal stress vector
            if this.fp
                figure(this.fig_strbar);
                quiver(this.x_gp,this.y_gp,this.fpx_gp,this.fpy_gp,'r');
                this.plotMesh(mdl,x,y,'k');
            end
            
            % Display flux results
            if this.fxx
                figure(this.fig_fxx);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.FXX_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.FXX_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.fyy
                figure(this.fig_fyy);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.FYY_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.FYY_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
            
            if this.fp
                figure(this.fig_fp);
                if this.smooth
                    this.plotNodeContourInplane(mdl,x,y,drv.Result.FP_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,x,y,drv.Result.FP_ELEMEXTRAP);
                end
                this.plotMesh(mdl,x,y,'k');
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Create figures for post-processing results.
        % Eight windows are created and positioned on the screen.
        % Each window plots a type of post-process result.
        function deform_fac = createFigsInplane(this,mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax)
            % Get current screen sizes.
            screen_sizes = get(0,'ScreenSize');
            
            % Create nodal displacement response data
            udispl = this.U(mdl.ID(1,:));
            vdispl = this.U(mdl.ID(2,:));
            udispl_min = min(udispl);
            udispl_max = max(udispl);
            vdispl_min = min(vdispl);
            vdispl_max = max(vdispl);
            
            % Compute deformed factor based on maximum displacement abs. value
            max_displ = max([abs(udispl_min), abs(udispl_max), abs(vdispl_min), abs(vdispl_max)]);
            min_plotsize = min(plot_xmax - plot_xmin,plot_ymax - plot_ymin);
            deform_fac = (min_plotsize / max_displ) * 0.15;
            
            % Create figure for mesh and deformed mesh plot and get its handle.
            % Locate figure at the up left corner of screen.
            this.fig_deform = figure;
            fig_deform_pos = get( this.fig_deform, 'Position' );
            fig_deform_pos(1) = 0;
            set( this.fig_deform, 'Position', fig_deform_pos );
            title_text = sprintf( 'Mesh and deformed mesh. Deformed factor: %s', num2str(deform_fac) );
            title( title_text );
            set( gca,'DataAspectRatio',[1 1 1] );
            shift = max_displ * deform_fac;
            plot_xmin_shift = plot_xmin - shift;
            plot_xmax_shift = plot_xmax + shift;
            plot_ymin_shift = plot_ymin - shift;
            plot_ymax_shift = plot_ymax + shift;
            axis([plot_xmin_shift plot_xmax_shift plot_ymin_shift plot_ymax_shift]);
            hold on;
            
            % Create figure for stress bar response plots and get handle to it.
            % Locate figure at the up right corner of screen.
            this.fig_strbar = figure;
            fig_strbar_pos = get( this.fig_strbar, 'Position' );
            fig_strbar_pos(1) = screen_sizes(3) - fig_strbar_pos(3);
            set( this.fig_strbar, 'Position', fig_strbar_pos );
            title( 'Principal stress directions' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            hold on;
            
            % Create figure for sigma x plot and get its handle.
            % Locate figure at the second level left side of screen.
            this.fig_sxx = figure;
            fig_sxx_pos = get( this.fig_sxx, 'Position' );
            fig_sxx_pos(1) = 0;
            fig_sxx_pos(2) = (screen_sizes(4) - fig_sxx_pos(4))/2;
            set( this.fig_sxx, 'Position', fig_sxx_pos );
            title( 'Sigma XX stress component' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.sxx_nodeextrap_min this.sxx_nodeextrap_max]);
            else
                caxis([this.sxx_elemextrap_min this.sxx_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for sigma y plot and get its handle.
            % Locate figure at the second level right side of screen.
            this.fig_syy = figure;
            fig_syy_pos = get( this.fig_syy, 'Position' );
            fig_syy_pos(1) = screen_sizes(3) - fig_syy_pos(3);
            fig_syy_pos(2) = (screen_sizes(4) - fig_syy_pos(4))/2;
            set( this.fig_syy, 'Position', fig_syy_pos );
            title( 'Sigma YY stress component' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.syy_nodeextrap_min this.syy_nodeextrap_max]);
            else
                caxis([this.syy_elemextrap_min this.syy_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for tau xy plot and get its handle.
            % Locate figure at the third level left side of screen.
            this.fig_txy = figure;
            fig_txy_pos = get( this.fig_txy, 'Position' );
            fig_txy_pos(1) = 0;
            fig_txy_pos(2) = (screen_sizes(4) - fig_txy_pos(4))/4;
            set( this.fig_txy, 'Position', fig_txy_pos );
            title( 'Tau XY stress component' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.txy_nodeextrap_min this.txy_nodeextrap_max]);
            else
                caxis([this.txy_elemextrap_min this.txy_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for sigma 1 plot and get its handle.
            % Locate figure at the third level right side of screen.
            this.fig_s1 = figure;
            fig_s1_pos = get( this.fig_s1, 'Position' );
            fig_s1_pos(1) = screen_sizes(3) - fig_s1_pos(3);
            fig_s1_pos(2) = (screen_sizes(4) - fig_s1_pos(4))/4;
            set( this.fig_s1, 'Position', fig_s1_pos );
            title( 'Maximum principal stress' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.s1_nodeextrap_min this.s1_nodeextrap_max]);
            else
                caxis([this.s1_elemextrap_min this.s1_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for sigma 2 plot and get its handle.
            % Locate figure at the forth level left side of screen.
            this.fig_s2 = figure;
            fig_s2_pos = get( this.fig_s2, 'Position' );
            fig_s2_pos(1) = 0;
            fig_s2_pos(2) = 0;
            set( this.fig_s2, 'Position', fig_s2_pos );
            title( 'Minimum principal stress' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.s2_nodeextrap_min this.s2_nodeextrap_max]);
            else
                caxis([this.s2_elemextrap_min this.s2_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for tau max. plot and get its handle.
            % Locate figure at the forth level right side of screen.
            this.fig_tmax = figure;
            fig_tmax_pos = get( this.fig_tmax, 'Position' );
            fig_tmax_pos(1) = screen_sizes(3) - fig_tmax_pos(3);
            fig_tmax_pos(2) = 0;
            set( this.fig_tmax, 'Position', fig_tmax_pos );
            title( 'Maximum shear stress' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.tmax_nodeextrap_min this.tmax_nodeextrap_max]);
            else
                caxis([this.tmax_elemextrap_min this.tmax_elemextrap_max]);
            end
            colorbar;
            hold on;
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Create figures for post-processing results.
        % Eight windows are created and positioned on the screen.
        % Each window plots a type of post-process result.
        function createFigsInplaneThermal(this,mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax)
            % Get current screen sizes.
            screen_sizes = get(0,'ScreenSize');
            
            % Create nodal displacement response data
            temperature = this.U(mdl.ID(1,:));
            temp_min = min(temperature);
            temp_max = max(temperature);
            
            % Create figure for temperature plot and get its handle.
            % Locate figure at the up left corner.
            this.fig_temp = figure;
            fig_temp_pos = get( this.fig_temp, 'Position' );
            fig_temp_pos(1) = 0;
            set( this.fig_temp, 'Position', fig_temp_pos );
            title( 'Temperature field' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            caxis([temp_min temp_max]);
            colorbar;
            hold on;
            
            % Create figure for sigma x plot and get its handle.
            % Locate figure at the second level left side of screen.
            this.fig_fxx = figure;
            fig_fxx_pos = get( this.fig_fxx, 'Position' );
            fig_fxx_pos(1) = 0;
            fig_fxx_pos(2) = (screen_sizes(4) - fig_fxx_pos(4))/2;
            set( this.fig_fxx, 'Position', fig_fxx_pos );
            title( 'Flux X component' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.fxx_nodeextrap_min this.fxx_nodeextrap_max]);
            else
                caxis([this.fxx_elemextrap_min this.fxx_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for sigma y plot and get its handle.
            % Locate figure at the third level left side of screen.
            this.fig_fyy = figure;
            fig_fyy_pos = get( this.fig_fyy, 'Position' );
            fig_fyy_pos(1) = 0;
            fig_fyy_pos(2) = (screen_sizes(4) - fig_fyy_pos(4))/4;
            set( this.fig_fyy, 'Position', fig_fyy_pos );
            title( 'Flux Y component' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.fyy_nodeextrap_min this.fyy_nodeextrap_max]);
            else
                caxis([this.fyy_elemextrap_min this.fyy_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for sigma 1 plot and get its handle.
            % Locate figure at the third level right side of screen.
            this.fig_fp = figure;
            fig_fp_pos = get( this.fig_fp, 'Position' );
            fig_fp_pos(1) = screen_sizes(3) - fig_fp_pos(3);
            fig_fp_pos(2) = (screen_sizes(4) - fig_fp_pos(4))/2;
            set( this.fig_fp, 'Position', fig_fp_pos );
            title( 'Flux in principal direction' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            if this.smooth
                caxis([this.fp_nodeextrap_min this.fp_nodeextrap_max]);
            else
                caxis([this.fp_elemextrap_min this.fp_elemextrap_max]);
            end
            colorbar;
            hold on;
            
            % Create figure for stress bar response plots and get handle to it.
            % Locate figure at the up right corner of screen.
            this.fig_strbar = figure;
            fig_strbar_pos = get( this.fig_strbar, 'Position' );
            fig_strbar_pos(1) = screen_sizes(3) - fig_strbar_pos(3);
            set( this.fig_strbar, 'Position', fig_strbar_pos );
            title( 'Flux principal directions' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            hold on;
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot mesh in current active figure.
        function plotMeshInplane(~,mdl,x,y,color)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            
            % Display mesh
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = x(node);
                    YY(j) = y(node);
                end
                node1 = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = x(node1);
                YY(nen+1) = y(node1);
                plot(XX,YY, color);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot deformed mesh in current active figure.
        function plotDeformMeshInplane(this,mdl,x,y,fct,color)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            
            u = this.U(mdl.ID(1,:));
            v = this.U(mdl.ID(2,:));
            
            % Display deformed mesh
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = x(node) + fct*u(node);
                    YY(j) = y(node) + fct*v(node);
                end
                node1 = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = x(node1) + fct*u(node1);
                YY(nen+1) = y(node1) + fct*v(node1);
                plot(XX,YY,color);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot current active element contour response in current active
        % figure.
        % Element contour is based on non average element results.
        function plotElemContourInplane(this,mdl,x,y,contour_type)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            % Display current active node contour response
            switch contour_type
                case drv.Result.SXX_ELEMEXTRAP
                    contour = this.sxx_elemextrap;
                case drv.Result.SYY_ELEMEXTRAP
                    contour = this.syy_elemextrap;
                case drv.Result.TXY_ELEMEXTRAP
                    contour = this.txy_elemextrap;
                case drv.Result.S1_ELEMEXTRAP
                    contour = this.s1_elemextrap;
                case drv.Result.S2_ELEMEXTRAP
                    contour = this.s2_elemextrap;
                case drv.Result.TMAX_ELEMEXTRAP
                    contour = this.tmax_elemextrap;
                case drv.Result.FXX_ELEMEXTRAP
                    contour = this.fxx_elemextrap;
                case drv.Result.FYY_ELEMEXTRAP
                    contour = this.fyy_elemextrap;
                case drv.Result.FP_ELEMEXTRAP
                    contour = this.fp_elemextrap;
            end
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    locid = mdl.elems(i).shape.ccwLocalNodeIds(j);
                    XX(j) = x(node);
                    YY(j) = y(node);
                    ZZ(j) = contour(locid,i);
                end
                node = mdl.elems(i).shape.ccwNodeIds(1);
                locid = mdl.elems(i).shape.ccwLocalNodeIds(1);
                XX(nen+1) = x(node);
                YY(nen+1) = y(node);
                ZZ(nen+1) = contour(locid,i);
                patch(XX,YY,ZZ);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot current active node contour response in current active
        % figure.
        % Node contour is based on nodal average adjacent element results.
        function plotNodeContourInplane(this,mdl,x,y,contour_type)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            % Display current active node contour response
            switch contour_type
                case drv.Result.SXX_NODEEXTRAP
                    contour = this.sxx_nodeextrap;
                case drv.Result.SYY_NODEEXTRAP
                    contour = this.syy_nodeextrap;
                case drv.Result.TXY_NODEEXTRAP
                    contour = this.txy_nodeextrap;
                case drv.Result.S1_NODEEXTRAP
                    contour = this.s1_nodeextrap;
                case drv.Result.S2_NODEEXTRAP
                    contour = this.s2_nodeextrap;
                case drv.Result.TMAX_NODEEXTRAP
                    contour = this.tmax_nodeextrap;
                case drv.Result.TEMP_NODE
                    contour = zeros(mdl.nnp,1);
                    for i = 1:mdl.nnp
                        contour(i) = this.U(mdl.ID(i));
                    end
                case drv.Result.FXX_NODEEXTRAP
                    contour = this.fxx_nodeextrap;
                case drv.Result.FYY_NODEEXTRAP
                    contour = this.fyy_nodeextrap;
                case drv.Result.FP_NODEEXTRAP
                    contour = this.fp_nodeextrap;
            end
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = x(node);
                    YY(j) = y(node);
                    ZZ(j) = contour(node);
                end
                node = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = x(node);
                YY(nen+1) = y(node);
                ZZ(nen+1) = contour(node);
                patch(XX,YY,ZZ);
                hold on;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Clear numerical garbage from analysis. (TEMPORARY SIMPLIFICATION !!)
        function clearSmallValues(this,mdl)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYM_STRESS
                this.clearSmallValuesInplane(mdl);
            elseif mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                   mdl.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.clearSmallValuesInplaneConduction(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Plot analysis results.
        function plot(this,mdl)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYM_STRESS
                this.plotInplane(mdl);
            elseif mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                   mdl.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.plotInplaneConduction(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Create figures for post-processing results.
        % Eight windows are created and positioned on the screen.
        % Each window plots a type of post-process result.
        function deform_fac = createFigs(this,mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYM_STRESS
                deform_fac = this.createFigsInplane(mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax);
            elseif mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                   mdl.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.createFigsInplaneThermal(mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax);
            end
        end
        
        %------------------------------------------------------------------
        % Plot mesh in current active figure.
        function plotMesh(this,mdl,x,y,color)
            if mdl.anm.type == fem.Anm.PLANE_STRESS     || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN     || ...
               mdl.anm.type == fem.Anm.AXISYM_STRESS    || ...
               mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
               mdl.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.plotMeshInplane(mdl,x,y,color);
            end
       end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in current active figure.
        function plotDeformMesh(this,mdl,x,y,fct,color)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYM_STRESS
                this.plotDeformMeshInplane(mdl,x,y,fct,color);
            end
        end
    end
end