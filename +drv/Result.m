%% Result Class
%
%% Description
%
% This class defines a result object in the StAnOOP program.
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
        smooth = logical(false);  % Smooth element results at common nodes
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
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        D = [];             % global displacement vector
        
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
        
        % Handles to plot figures
        fig_deform = [];    % figure for mesh and deformed mesh plot
        fig_strbar = [];    % figure for stress bar response plots
        fig_sxx    = [];    % figure for sigma x plot
        fig_syy    = [];    % figure for sigma y plot
        fig_txy    = [];    % figure for tau xy plot
        fig_s1     = []     % figure for sigma 1 plot
        fig_s2     = []     % figure for sigma 2 plot
        fig_tmax   = []     % figure for tau max. plot
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function res = Result()
            return;
        end
    end
    
    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Clear numerical garbage from analysis. (TEMPORARY SIMPLIFICATION !!)
        function clearSmallValuesInplane(res,mdl)
            if (abs(res.sxx_gp_min) < 0.00001 && abs(res.sxx_gp_max) < 0.00001)
                res.sxx_gp_min = 0.0;
                res.sxx_gp_max = 0.0;
                res.sxx_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(res.syy_gp_min) < 0.00001 && abs(res.syy_gp_max) < 0.00001)
                res.syy_gp_min = 0.0;
                res.syy_gp_max = 0.0;
                res.syy_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(res.txy_gp_min) < 0.00001 && abs(res.txy_gp_max) < 0.00001)
                res.txy_gp_min = 0.0;
                res.txy_gp_max = 0.0;
                res.txy_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(res.s1_gp_min) < 0.00001 && abs(res.s1_gp_max) < 0.00001)
                res.s1_gp_min = 0.0;
                res.s1_gp_max = 0.0;
                res.s1_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(res.s2_gp_min) < 0.00001 && abs(res.s2_gp_max) < 0.00001)
                res.s2_gp_min = 0.0;
                res.s2_gp_max = 0.0;
                res.s2_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(res.tmax_gp_min) < 0.00001 && abs(res.tmax_gp_max) < 0.00001)
                res.tmax_gp_min = 0.0;
                res.tmax_gp_max = 0.0;
                res.tmax_gp = zeros(mdl.elems(1).gstress_npts,mdl.nel);
            end
            if (abs(res.sxx_elemextrap_min) < 0.00001 && abs(res.sxx_elemextrap_max) < 0.00001)
                res.sxx_elemextrap_min = 0.0;
                res.sxx_elemextrap_max = 0.0;
                res.sxx_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(res.syy_elemextrap_min) < 0.00001 && abs(res.syy_elemextrap_max) < 0.00001)
                res.syy_elemextrap_min = 0.0;
                res.syy_elemextrap_max = 0.0;
                res.syy_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(res.txy_elemextrap_min) < 0.00001 && abs(res.txy_elemextrap_max) < 0.00001)
                res.txy_elemextrap_min = 0.0;
                res.txy_elemextrap_max = 0.0;
                res.txy_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(res.s1_elemextrap_min) < 0.00001 && abs(res.s1_elemextrap_max) < 0.00001)
                res.s1_elemextrap_min = 0.0;
                res.s1_elemextrap_max = 0.0;
                res.s1_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(res.s2_elemextrap_min) < 0.00001 && abs(res.s2_elemextrap_max) < 0.00001)
                res.s2_elemextrap_min = 0.0;
                res.s2_elemextrap_max = 0.0;
                res.s2_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
            if (abs(res.tmax_elemextrap_min) < 0.00001 && abs(res.tmax_elemextrap_max) < 0.00001)
                res.tmax_elemextrap_min = 0.0;
                res.tmax_elemextrap_max = 0.0;
                res.tmax_elemextrap = zeros(mdl.elems(1).shape.nen,mdl.nel);
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Plot analysis results.
        function plotInplane(res,mdl)
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
            
            % Create figures (windows) for displaying results.
            deform_fac = res.createFigs(mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax);
            
            % Display deformed mesh.
            figure(res.fig_deform);
            res.plotMesh(mdl,x,y,'k');
            res.plotDeformMesh(mdl,x,y,deform_fac,'b');
            
            % Display principal stress vectors.
            figure(res.fig_strbar);
            res.plotMesh(mdl,x,y,'k');
            quiver(res.x_gp,res.y_gp,res.s1x_gp,res.s1y_gp,'r');
            quiver(res.x_gp,res.y_gp,res.s2x_gp,res.s2y_gp,'b');
            
            % Display stress results.
            if res.sxx
                figure(res.fig_sxx);
                if res.smooth
                    res.plotNodeContourInplane(mdl,x,y,drv.Result.SXX_NODEEXTRAP);
                else
                    res.plotElemContourInplane(mdl,x,y,drv.Result.SXX_ELEMEXTRAP);
                end
                res.plotMesh(mdl,x,y,'k');
            end
            
            if res.syy
                figure(res.fig_syy);
                if res.smooth
                    res.plotNodeContourInplane(mdl,x,y,drv.Result.SYY_NODEEXTRAP);
                else
                    res.plotElemContourInplane(mdl,x,y,drv.Result.SYY_ELEMEXTRAP);
                end
                res.plotMesh(mdl,x,y,'k');
            end
            
            if res.txy
                figure(res.fig_txy);
                if res.smooth
                    res.plotNodeContourInplane(mdl,x,y,drv.Result.TXY_NODEEXTRAP);
                else
                    res.plotElemContourInplane(mdl,x,y,drv.Result.TXY_ELEMEXTRAP);
                end
                res.plotMesh(mdl,x,y,'k');
            end
            
            if res.s1
                figure(res.fig_s1);
                if res.smooth
                    res.plotNodeContourInplane(mdl,x,y,drv.Result.S1_NODEEXTRAP);
                else
                    res.plotElemContourInplane(mdl,x,y,drv.Result.S1_ELEMEXTRAP);
                end
                res.plotMesh(mdl,x,y,'k');
            end
            
            if res.s2
                figure(res.fig_s2);
                if res.smooth
                    res.plotNodeContourInplane(mdl,x,y,drv.Result.S2_NODEEXTRAP);
                else
                    res.plotElemContourInplane(mdl,x,y,drv.Result.S2_ELEMEXTRAP);
                end
                res.plotMesh(mdl,x,y,'k');
            end
            
            if res.taumax
                figure(res.fig_tmax);
                if res.smooth
                    res.plotNodeContourInplane(mdl,x,y,drv.Result.TMAX_NODEEXTRAP);
                else
                    res.plotElemContourInplane(mdl,x,y,drv.Result.TMAX_ELEMEXTRAP);
                end
                res.plotMesh(mdl,x,y,'k');
            end
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Create figures for post-processing results.
        % Eight windows are created and positioned on the screen.
        % Each window plots a type of post-process result.
        function deform_fac = createFigsInplane(res,mdl, ...
                                   plot_xmin,plot_xmax,plot_ymin,plot_ymax)
            % Get current screen sizes.
            screen_sizes = get(0,'ScreenSize');
            
            % Create nodal displacement response data
            udispl = res.D(mdl.ID(1,:));
            vdispl = res.D(mdl.ID(2,:));
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
            res.fig_deform = figure;
            fig_deform_pos = get( res.fig_deform, 'Position' );
            fig_deform_pos(1) = 0;
            set( res.fig_deform, 'Position', fig_deform_pos );
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
            res.fig_strbar = figure;
            fig_strbar_pos = get( res.fig_strbar, 'Position' );
            fig_strbar_pos(1) = screen_sizes(3) - fig_strbar_pos(3);
            set( res.fig_strbar, 'Position', fig_strbar_pos );
            title( 'Principal stress directions' );
            set( gca,'DataAspectRatio',[1 1 1] );
            axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
            hold on;
            
             % Create figure for sigma x plot and get its handle.
             % Locate figure at the second level left side of screen.
             res.fig_sxx = figure;
             fig_sxx_pos = get( res.fig_sxx, 'Position' );
             fig_sxx_pos(1) = 0;
             fig_sxx_pos(2) = (screen_sizes(4) - fig_sxx_pos(4))/2;
             set( res.fig_sxx, 'Position', fig_sxx_pos );
             title( 'Sigma XX stress component' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             if res.smooth
                 caxis([res.sxx_nodeextrap_min res.sxx_nodeextrap_max]);
             else
                 caxis([res.sxx_elemextrap_min res.sxx_elemextrap_max]);
             end
             colorbar;
             hold on;
            
             % Create figure for sigma y plot and get its handle.
             % Locate figure at the second level right side of screen.
             res.fig_syy = figure;
             fig_syy_pos = get( res.fig_syy, 'Position' );
             fig_syy_pos(1) = screen_sizes(3) - fig_syy_pos(3);
             fig_syy_pos(2) = (screen_sizes(4) - fig_syy_pos(4))/2;
             set( res.fig_syy, 'Position', fig_syy_pos );
             title( 'Sigma YY stress component' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             if res.smooth
                 caxis([res.syy_nodeextrap_min res.syy_nodeextrap_max]);
             else
                 caxis([res.syy_elemextrap_min res.syy_elemextrap_max]);
             end
             colorbar;
             hold on;
             
             % Create figure for tau xy plot and get its handle.
             % Locate figure at the third level left side of screen.
             res.fig_txy = figure;
             fig_txy_pos = get( res.fig_txy, 'Position' );
             fig_txy_pos(1) = 0;
             fig_txy_pos(2) = (screen_sizes(4) - fig_txy_pos(4))/4;
             set( res.fig_txy, 'Position', fig_txy_pos );
             title( 'Tau XY stress component' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             if res.smooth
                 caxis([res.txy_nodeextrap_min res.txy_nodeextrap_max]);
             else
                 caxis([res.txy_elemextrap_min res.txy_elemextrap_max]);
             end
             colorbar;
             hold on;
             
             % Create figure for sigma 1 plot and get its handle.
             % Locate figure at the third level right side of screen.
             res.fig_s1 = figure;
             fig_s1_pos = get( res.fig_s1, 'Position' );
             fig_s1_pos(1) = screen_sizes(3) - fig_s1_pos(3);
             fig_s1_pos(2) = (screen_sizes(4) - fig_s1_pos(4))/4;
             set( res.fig_s1, 'Position', fig_s1_pos );
             title( 'Maximum principal stress' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             if res.smooth
                 caxis([res.s1_nodeextrap_min res.s1_nodeextrap_max]);
             else
                 caxis([res.s1_elemextrap_min res.s1_elemextrap_max]);
             end
             colorbar;
             hold on;
             
             % Create figure for sigma 2 plot and get its handle.
             % Locate figure at the forth level left side of screen.
             res.fig_s2 = figure;
             fig_s2_pos = get( res.fig_s2, 'Position' );
             fig_s2_pos(1) = 0;
             fig_s2_pos(2) = 0;
             set( res.fig_s2, 'Position', fig_s2_pos );
             title( 'Minimum principal stress' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             if res.smooth
                 caxis([res.s2_nodeextrap_min res.s2_nodeextrap_max]);
             else
                 caxis([res.s2_elemextrap_min res.s2_elemextrap_max]);
             end
             colorbar;
             hold on;
             
             % Create figure for tau max. plot and get its handle.
             % Locate figure at the forth level right side of screen.
             res.fig_tmax = figure;
             fig_tmax_pos = get( res.fig_tmax, 'Position' );
             fig_tmax_pos(1) = screen_sizes(3) - fig_tmax_pos(3);
             fig_tmax_pos(2) = 0;
             set( res.fig_tmax, 'Position', fig_tmax_pos );
             title( 'Maximum shear stress' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             if res.smooth
                 caxis([res.tmax_nodeextrap_min res.tmax_nodeextrap_max]);
             else
                 caxis([res.tmax_elemextrap_min res.tmax_elemextrap_max]);
             end
             colorbar;
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
        function plotDeformMeshInplane(res,mdl,x,y,fct,color)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            
            u = res.D(mdl.ID(1,:));
            v = res.D(mdl.ID(2,:));
            
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
        function plotElemContourInplane(res,mdl,x,y,contour_type)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            % Display current active node contour response
            switch contour_type
                case drv.Result.SXX_ELEMEXTRAP
                    contour = res.sxx_elemextrap;
                case drv.Result.SYY_ELEMEXTRAP
                    contour = res.syy_elemextrap;
                case drv.Result.TXY_ELEMEXTRAP
                    contour = res.txy_elemextrap;
                case drv.Result.S1_ELEMEXTRAP
                    contour = res.s1_elemextrap;
                case drv.Result.S2_ELEMEXTRAP
                    contour = res.s2_elemextrap;
                case drv.Result.TMAX_ELEMEXTRAP
                    contour = res.tmax_elemextrap;
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
        function plotNodeContourInplane(res,mdl,x,y,contour_type)
            maxNen = mdl.maxNumElemNodes;
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            % Display current active node contour response
            switch contour_type
                case drv.Result.SXX_NODEEXTRAP
                    contour = res.sxx_nodeextrap;
                case drv.Result.SYY_NODEEXTRAP
                    contour = res.syy_nodeextrap;
                case drv.Result.TXY_NODEEXTRAP
                    contour = res.txy_nodeextrap;
                case drv.Result.S1_NODEEXTRAP
                    contour = res.s1_nodeextrap;
                case drv.Result.S2_NODEEXTRAP
                    contour = res.s2_nodeextrap;
                case drv.Result.TMAX_NODEEXTRAP
                    contour = res.tmax_nodeextrap;
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
        function clearSmallValues(res,mdl)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                res.clearSmallValuesInplane(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Plot analysis results.
        function plot(res,mdl)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                res.plotInplane(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Create figures for post-processing results.
        % Eight windows are created and positioned on the screen.
        % Each window plots a type of post-process result.
        function deform_fac = createFigs(res,mdl, ...
                                   plot_xmin,plot_xmax,plot_ymin,plot_ymax)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                deform_fac = res.createFigsInplane(mdl, ...
                                  plot_xmin,plot_xmax,plot_ymin,plot_ymax);
            end
        end
        
        %------------------------------------------------------------------
        % Plot mesh in current active figure.
        function plotMesh(res,mdl,x,y,color)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                res.plotMeshInplane(mdl,x,y,color);
            end
       end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in current active figure.
        function plotDeformMesh(res,mdl,x,y,fct,color)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                res.plotDeformMeshInplane(mdl,x,y,fct,color);
            end
        end
    end
end