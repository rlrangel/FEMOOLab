%% Plot Class
%
%% Description
%
% This class defines a plotting object in the FEMOOLab program.
% A plotting object is responsible for displaying the analysis results in
% Matlab figure windows.
%
% This class currently deals with the plotting of all types of analysis
% models. However, for the sake of oraganization, it should become a super-
% class with a plotting sub-class for each analysis model.
%
classdef Plot < handle
    %% Constant values for contour types
    properties (Constant = true, Access = public)
        % Structural analysis node results
        DX_NODE           = int32(1);
        DY_NODE           = int32(2);
        DZ_NODE           = int32(3);
        RX_NODE           = int32(4);
        RY_NODE           = int32(5);
        RZ_NODE           = int32(6);
        
        % Structural analysis element extrapolated results
        SXX_ELEMEXTRAP    = int32(7);
        SYY_ELEMEXTRAP    = int32(8);
        SZZ_ELEMEXTRAP    = int32(9);
        TXY_ELEMEXTRAP    = int32(10);
        TXZ_ELEMEXTRAP    = int32(11);
        TYZ_ELEMEXTRAP    = int32(12);
        S1_ELEMEXTRAP     = int32(13);
        S2_ELEMEXTRAP     = int32(14);
        S3_ELEMEXTRAP     = int32(15);
        TMAX_ELEMEXTRAP   = int32(16);
        QXZ_ELEMEXTRAP    = int32(17);
        QYZ_ELEMEXTRAP    = int32(18);
        MXX_ELEMEXTRAP    = int32(19);
        MYY_ELEMEXTRAP    = int32(20);
        MXY_ELEMEXTRAP    = int32(21);
        M1_ELEMEXTRAP     = int32(22);
        M2_ELEMEXTRAP     = int32(23);
        TORMAX_ELEMEXTRAP = int32(24);
        
        % Structural analysis node smoothed results
        SXX_NODEEXTRAP    = int32(25);
        SYY_NODEEXTRAP    = int32(26);
        SZZ_NODEEXTRAP    = int32(27);
        TXY_NODEEXTRAP    = int32(28);
        TXZ_NODEEXTRAP    = int32(29);
        TYZ_NODEEXTRAP    = int32(30);
        S1_NODEEXTRAP     = int32(31);
        S2_NODEEXTRAP     = int32(32);
        S3_NODEEXTRAP     = int32(33);
        TMAX_NODEEXTRAP   = int32(34);
        QXZ_NODEEXTRAP    = int32(35);
        QYZ_NODEEXTRAP    = int32(36);
        MXX_NODEEXTRAP    = int32(37);
        MYY_NODEEXTRAP    = int32(38);
        MXY_NODEEXTRAP    = int32(39);
        M1_NODEEXTRAP     = int32(40);
        M2_NODEEXTRAP     = int32(41);
        TORMAX_NODEEXTRAP = int32(42);
        
        % Thermal analysis node results
        TEMP_NODE         = int32(43);
        
        % Thermal analysis element extrapolated results
        FXX_ELEMEXTRAP    = int32(44);
        FYY_ELEMEXTRAP    = int32(45);
        FZZ_ELEMEXTRAP    = int32(46);
        FM_ELEMEXTRAP     = int32(47);
        
        % Thermal analysis node smoothed results
        FXX_NODEEXTRAP    = int32(48);
        FYY_NODEEXTRAP    = int32(49);
        FZZ_NODEEXTRAP    = int32(50);
        FM_NODEEXTRAP     = int32(51);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Handles to plot figures
        fig_lbl    = [];  % figure for mesh labels
        fig_deform = [];  % figure for mesh and deformed mesh plot
        fig_dx     = [];  % figure for displacement x plot
        fig_dy     = [];  % figure for displacement y plot
        fig_dz     = [];  % figure for displacement z plot
        fig_rx     = [];  % figure for rotation x plot
        fig_ry     = [];  % figure for rotation y plot
        fig_sxx    = [];  % figure for sigma x plot
        fig_syy    = [];  % figure for sigma y plot
        fig_txy    = [];  % figure for tau xy plot
        fig_s1     = [];  % figure for sigma 1 plot
        fig_s2     = [];  % figure for sigma 2 plot
        fig_tmax   = [];  % figure for tau max. plot
        fig_prcstr = [];  % figure for principal stress directions
        fig_qxz    = [];  % figure for shear xz plot
        fig_qyz    = [];  % figure for shear yz plot
        fig_mxx    = [];  % figure for moment x plot
        fig_myy    = [];  % figure for moment y plot
        fig_mxy    = [];  % figure for moment xy plot
        fig_m1     = [];  % figure for moment 1 plot
        fig_m2     = [];  % figure for moment 2 plot
        fig_tormax = [];  % figure for torsion max plot
        fig_prcmom = [];  % figure for principal moment directions
        fig_temp   = [];  % figure for temperature field
        fig_fxx    = [];  % figure for flux x plot
        fig_fyy    = [];  % figure for flux x plot
        fig_fm     = [];  % figure for flux module plot
        fig_prcflx = [];  % figure for flux directions
        
        % Bounding box for plot figures
        x_coord    double = double.empty;
        y_coord    double = double.empty;
        plot_xmin  double = double.empty;
        plot_xmax  double = double.empty;
        plot_ymin  double = double.empty;
        plot_ymax  double = double.empty;
        
        % Visual properties
        deform_fac   double = double.empty;   % Scale factor for deformed mesh
        color_mesh   = 'k';                   % mesh color
        color_deform = 'b';                   % deformed mesh color
        fcn_dataTip  = @drv.cb_dataTipCursor; % handle to function to manage data tip showing
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Plot(sim)
            if (nargin > 0)
                this.setupBoundingBox(sim.mdl);
                
                if (sim.anl.type == fem.Anl.LINEAR_STATIC)
                    this.plotStatic(sim.mdl);
                elseif (sim.anl.type == fem.Anl.LINEAR_TRANSIENT)
                    this.plotTransientCurves(sim.mdl);
                    this.plotTransientContours(sim.mdl,sim.anl);
                end
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Setup bounding box for plot figures.
        function setupBoundingBox(this,mdl)
            % Assemble vectors of nodal coordinates
            x = zeros(mdl.nnp,1);
            y = zeros(mdl.nnp,1);
            for i = 1:mdl.nnp
                x(i) = mdl.nodes(i).coord(1);
                y(i) = mdl.nodes(i).coord(2);
            end
            this.x_coord = x;
            this.y_coord = y;
            
            % Setup bounding box for displaying results
            min_x = min(x);
            max_x = max(x);
            min_y = min(y);
            max_y = max(y);
            size_x = max_x-min_x;
            size_y = max_y-min_y;
            cx = min_x + size_x * 0.5;
            cy = min_y + size_y * 0.5;
            this.plot_xmin = cx - size_x * 0.55;
            this.plot_xmax = cx + size_x * 0.55;
            this.plot_ymin = cy - size_y * 0.55;
            this.plot_ymax = cy + size_y * 0.55;
        end
        
        %------------------------------------------------------------------
        % Plot static analysis results.
        function plotStatic(this,mdl)
            % Create figures (windows) for displaying results
            this.createStaticFigs(mdl);
            
            % Display results
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
                this.plotMesh(mdl);
                this.plotMeshLabels(mdl);
            end
            
            if (mdl.res.deform)
                figure(this.fig_deform);
                this.plotMesh(mdl);
                this.plotDeformMesh(mdl);
            end
            
            if (mdl.res.dx)
                figure(this.fig_dx);
                this.plotNodeContour(mdl,this.DX_NODE);
            end
            
            if (mdl.res.dy)
                figure(this.fig_dy);
                this.plotNodeContour(mdl,this.DY_NODE);
            end
            
            if (mdl.res.dz)
                figure(this.fig_dz);
                this.plotNodeContour(mdl,this.DZ_NODE);
            end
            
            if (mdl.res.rx)
                figure(this.fig_rx);
                this.plotNodeContour(mdl,this.RX_NODE);
            end
            
            if (mdl.res.ry)
                figure(this.fig_ry);
                this.plotNodeContour(mdl,this.RY_NODE);
            end
            
            if (mdl.res.sxx)
                figure(this.fig_sxx);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.SXX_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.SXX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.syy)
                figure(this.fig_syy);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.SYY_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.SYY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.txy)
                figure(this.fig_txy);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.TXY_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.TXY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.s1)
                figure(this.fig_s1);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.S1_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.S1_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.s2)
                figure(this.fig_s2);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.S2_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.S2_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.taumax)
                figure(this.fig_tmax);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.TMAX_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.TMAX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.s1 || mdl.res.s2)
                figure(this.fig_prcstr);
                this.plotMesh(mdl);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.s1x_gp,mdl.res.s1y_gp,'r');
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.s2x_gp,mdl.res.s2y_gp,'b');
            end
            
            if (mdl.res.qxz)
                figure(this.fig_qxz);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.QXZ_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.QXZ_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.qyz)
                figure(this.fig_qyz);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.QYZ_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.QYZ_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.mxx)
                figure(this.fig_mxx);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.MXX_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.MXX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.myy)
                figure(this.fig_myy);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.MYY_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.MYY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.mxy)
                figure(this.fig_mxy);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.MXY_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.MXY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.m1)
                figure(this.fig_m1);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.M1_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.M1_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.m2)
                figure(this.fig_m2);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.M2_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.M2_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.tormax)
                figure(this.fig_tormax);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.TORMAX_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.TORMAX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.m1 || mdl.res.m2)
                figure(this.fig_prcmom);
                this.plotMesh(mdl);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.m1x_gp,mdl.res.m1y_gp,'r');
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.m2x_gp,mdl.res.m2y_gp,'b');
            end
            
            if (mdl.res.temp)
                figure(this.fig_temp);
                this.plotNodeContour(mdl,this.TEMP_NODE);
            end
            
            if (mdl.res.fxx)
                figure(this.fig_fxx);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.FXX_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.FXX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.fyy)
                figure(this.fig_fyy);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.FYY_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.FYY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.fm)
                figure(this.fig_fm);
                if (mdl.res.smooth)
                    this.plotNodeContour(mdl,this.FM_NODEEXTRAP);
                else
                    this.plotElemContour(mdl,this.FM_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.fm)
                figure(this.fig_prcflx);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.fmx_gp,mdl.res.fmy_gp,'r');
                this.plotMesh(mdl);
            end
            
            % Set order to display (reverse)
            if (mdl.res.tormax)
                figure(this.fig_tormax);
            end
            if (mdl.res.m2)
                figure(this.fig_m2);
            end
            if (mdl.res.m1)
                figure(this.fig_m1);
            end
            if (mdl.res.m1 || mdl.res.m2)
                figure(this.fig_prcmom);
            end
            if (mdl.res.mxy)
                figure(this.fig_mxy);
            end
            if (mdl.res.myy)
                figure(this.fig_myy);
            end
            if (mdl.res.mxx)
                figure(this.fig_mxx);
            end
            if (mdl.res.qyz)
                figure(this.fig_qyz);
            end
            if (mdl.res.qxz)
                figure(this.fig_qxz);
            end
            if (mdl.res.taumax)
                figure(this.fig_tmax);
            end
            if (mdl.res.s2)
                figure(this.fig_s2);
            end
            if (mdl.res.s1)
                figure(this.fig_s1);
            end
            if (mdl.res.s1 || mdl.res.s2)
                figure(this.fig_prcstr);
            end
            if (mdl.res.txy)
                figure(this.fig_txy);
            end
            if (mdl.res.syy)
                figure(this.fig_syy);
            end
            if (mdl.res.sxx)
                figure(this.fig_sxx);
            end
            if (mdl.res.ry)
                figure(this.fig_ry);
            end
            if (mdl.res.rx)
                figure(this.fig_rx);
            end
            if (mdl.res.dz)
                figure(this.fig_dz);
            end
            if (mdl.res.dy)
                figure(this.fig_dy);
            end
            if (mdl.res.dx)
                figure(this.fig_dx);
            end
            if (mdl.res.deform)
                figure(this.fig_deform);
            end
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
            end
            
            if (mdl.res.fyy)
                figure(this.fig_fyy);
            end
            if (mdl.res.fxx)
                figure(this.fig_fxx);
            end
            if (mdl.res.fm)
                figure(this.fig_fm);
                figure(this.fig_prcflx);
            end
            if (mdl.res.temp)
                figure(this.fig_temp);
            end
        end
        
        %------------------------------------------------------------------
        % Plot transient analysis contours.
        function plotTransientContours(this,mdl,anl)
            % Display mesh labels
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                this.fig_lbl = figure;
                ax = gca;
                movegui(ax,'center')
                set(ax,'DataAspectRatio',[1 1 1]);
                title('Mesh labels');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
                this.plotMesh(mdl);
                this.plotMeshLabels(mdl);
            end
            
            % Pot transient animation
            if (mdl.res.temp)
                % Number of frames
                steps = 1:mdl.res.steps;
                steps = steps(rem(steps,mdl.res.output_freq)==0);
                loops = length(steps);
                
                % Create arrays with data to be ploted in all steps
                if (isempty(mdl.res.maxNen))
                    maxNen = mdl.maxNumElemNodes();
                else
                    maxNen = mdl.res.maxNen;
                end
                XX = zeros(mdl.nel,maxNen+1,loops);
                YY = zeros(mdl.nel,maxNen+1,loops);
                ZZ = zeros(mdl.nel,maxNen+1,loops);
                contour = mdl.res.U(mdl.ID(1,:),steps);
                for i = 1:mdl.nel
                    nen = mdl.elems(i).shape.nen;
                    for j = 1:nen
                        node= mdl.elems(i).shape.ccwNodeIds(j);
                        XX(i,j,:) = this.x_coord(node);
                        YY(i,j,:) = this.y_coord(node);
                        ZZ(i,j,:) = contour(node,:);
                    end
                    node = mdl.elems(i).shape.ccwNodeIds(1);
                    XX(i,nen+1,:) = this.x_coord(node);
                    YY(i,nen+1,:) = this.y_coord(node);
                    ZZ(i,nen+1,:) = contour(node,:);
                end
                
                % Get max and min temperature values
                temp_min = min(min(mdl.res.U));
                temp_max = max(max(mdl.res.U));
                
                % Create figure object
                this.fig_temp = figure;
                ax = gca;
                movegui(ax,'center');
                set(ax,'DataAspectRatio',[1 1 1],'Colormap',jet);
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (temp_min < temp_max)
                    caxis([temp_min temp_max]);
                end
                colorbar;
                
                % Plot animation
                for i = 1:loops
                    for j = 1:mdl.nel
                        patch(XX(j,:,i),YY(j,:,i),ZZ(j,:,i));
                    end
                    title(['Temperature, t = ',num2str(anl.incr*steps(i))],'FontSize',12)
                    pause(eps);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plot transient analysis curves.
        function plotTransientCurves(~,mdl)
            
            % Nodal temperature
            if (~isempty(mdl.res.curve_temp))
                ncurves = length(mdl.res.curve_temp);
                leg = strings(1,ncurves);
                f1 = figure; hold on; grid on;
                f2 = figure; hold on; grid on;
                for i = 1:ncurves
                    dof = mdl.ID(1,mdl.res.curve_temp(i));
                    x  = mdl.res.times;
                    y  = mdl.res.U(dof,:);
                    yt = mdl.res.Ut(dof,:);
                    legend_txt = sprintf('Node %d',mdl.res.curve_temp(i));
                    leg(i) = legend_txt;
                    % Plot temperature values
                    figure(f1);
                    plot(x,y);
                    % Plot temperature rate of cahnge
                    figure(f2);
                    plot(x,yt);
                end
                figure(f2);
                xlabel('Time');
                ylabel('Temperature Rate of Change');
                title('Nodal Temperature Rate of Change');
                legend(leg);
                movegui(gca,'center')
                
                figure(f1);
                xlabel('Time');
                ylabel('Temperature');
                title('Nodal Temperature');
                legend(leg);
                movegui(gca,'center')
            end
            
        end
        
        %------------------------------------------------------------------
        % Create figures for plotting 2D inplane static analysis results.
        function createStaticFigsInplane(this,mdl)
            % Mesh labels
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                this.fig_lbl = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Mesh labels');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
            
            % Deformed configuration
            if (mdl.res.deform)
                % Create nodal displacement response data
                udispl     = mdl.res.U(mdl.ID(1,:));
                vdispl     = mdl.res.U(mdl.ID(2,:));
                udispl_min = min(udispl);
                udispl_max = max(udispl);
                vdispl_min = min(vdispl);
                vdispl_max = max(vdispl);
                
                % Compute deformed factor based on maximum displacement abs. value
                max_displ = max([abs(udispl_min), abs(udispl_max), abs(vdispl_min), abs(vdispl_max)]);
                min_plotsize = min(this.plot_xmax - this.plot_xmin, this.plot_ymax - this.plot_ymin);
                
                if (mdl.res.scl == 0)
                    this.deform_fac = (min_plotsize / max_displ) * 0.15;
                else
                    this.deform_fac = mdl.res.scl;
                end
                
                % Create figure for mesh and deformed mesh plot
                this.fig_deform = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title_text = sprintf('Mesh and deformed mesh. Deformed factor: %s',num2str(this.deform_fac));
                title(title_text);
                shift = max_displ * this.deform_fac;
                plot_xmin_shift = this.plot_xmin - shift;
                plot_xmax_shift = this.plot_xmax + shift;
                plot_ymin_shift = this.plot_ymin - shift;
                plot_ymax_shift = this.plot_ymax + shift;
                axis([plot_xmin_shift plot_xmax_shift plot_ymin_shift plot_ymax_shift]);
                hold on;
            end
            
            % Displacement X
            if (mdl.res.dx)
                % Create figure for displacement plot
                this.fig_dx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Displacement X');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Displacement Y
            if (mdl.res.dy)
                % Create figure for displacement plot
                this.fig_dy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Displacement Y');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Sigma X
            if (mdl.res.sxx)
                this.fig_sxx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Sigma XX stress component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.sxx_nodeextrap_min mdl.res.sxx_nodeextrap_max]);
                else
                    caxis([mdl.res.sxx_elemextrap_min mdl.res.sxx_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Sigma Y
            if (mdl.res.syy)
                this.fig_syy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Sigma YY stress component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.syy_nodeextrap_min mdl.res.syy_nodeextrap_max]);
                else
                    caxis([mdl.res.syy_elemextrap_min mdl.res.syy_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Tau XY
            if (mdl.res.txy)
                this.fig_txy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Tau XY stress component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.txy_nodeextrap_min mdl.res.txy_nodeextrap_max]);
                else
                    caxis([mdl.res.txy_elemextrap_min mdl.res.txy_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % sigma 1
            if (mdl.res.s1)
                this.fig_s1 = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum principal stress');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.s1_nodeextrap_min mdl.res.s1_nodeextrap_max]);
                else
                    caxis([mdl.res.s1_elemextrap_min mdl.res.s1_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % sigma 2
            if (mdl.res.s2)
                this.fig_s2 = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Minimum principal stress');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.s2_nodeextrap_min mdl.res.s2_nodeextrap_max]);
                else
                    caxis([mdl.res.s2_elemextrap_min mdl.res.s2_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Tau max
            if (mdl.res.taumax)
                this.fig_tmax = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum shear stress');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.tmax_nodeextrap_min mdl.res.tmax_nodeextrap_max]);
                else
                    caxis([mdl.res.tmax_elemextrap_min mdl.res.tmax_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Principal stress directions
            if (mdl.res.s1 || mdl.res.s2)
                this.fig_prcstr = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Principal stress directions');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
            
            % Temperature field
            if (mdl.res.temp)
                % Create figure for temperature plot
                this.fig_temp = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Temperature field');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Heat flux X
            if (mdl.res.fxx)
                this.fig_fxx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Heat flux X component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.fxx_nodeextrap_min mdl.res.fxx_nodeextrap_max]);
                else
                    caxis([mdl.res.fxx_elemextrap_min mdl.res.fxx_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Heat flux Y
            if (mdl.res.fyy)
                this.fig_fyy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Heat flux Y component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.fyy_nodeextrap_min mdl.res.fyy_nodeextrap_max]);
                else
                    caxis([mdl.res.fyy_elemextrap_min mdl.res.fyy_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Heat flux module
            if (mdl.res.fm)
                this.fig_fm = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Heat flux module');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.fm_nodeextrap_min mdl.res.fm_nodeextrap_max]);
                else
                    caxis([mdl.res.fm_elemextrap_min mdl.res.fm_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Heat flux directions
            if (mdl.res.fm)
                this.fig_prcflx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Heat flux directions');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Create figures for plotting 2D outplane static analysis results.
        function createStaticFigsOutplane(this,mdl)
            % Mesh labels
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                this.fig_lbl = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Mesh labels');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
            
            % Deformed configuration
            if (mdl.res.deform)
                % Create nodal displacement response data
                wdispl     = mdl.res.U(mdl.ID(1,:));
                wdispl_min = min(wdispl);
                wdispl_max = max(wdispl);
                
                % Compute deformed factor based on maximum displacement abs. value
                max_displ = max([abs(wdispl_min), abs(wdispl_max)]);
                min_plotsize = min(this.plot_xmax - this.plot_xmin, this.plot_ymax - this.plot_ymin);
                
                if (mdl.res.scl == 0)
                    this.deform_fac = (min_plotsize / max_displ) * 0.15;
                else
                    this.deform_fac = mdl.res.scl;
                end
                
                % Create figure for mesh and deformed mesh plot
                this.fig_deform = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title_text = sprintf('Mesh and deformed mesh. Deformed factor: %s',num2str(this.deform_fac));
                title(title_text);
                shift = max_displ * this.deform_fac;
                plot_zmin = -shift;
                plot_zmax =  shift;
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax plot_zmin plot_zmax]);
                view(3);
                hold on;
            end
            
            % Displacement Z
            if (mdl.res.dz)
                % Create figure for displacement plot
                this.fig_dz = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Displacement Z');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Rotation X
            if (mdl.res.rx)
                % Create figure for rotation plot
                this.fig_rx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Rotation X');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Rotation Y
            if (mdl.res.ry)
                % Create figure for rotation plot
                this.fig_ry = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Rotation Y');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Shear XZ
            if (mdl.res.qxz)
                this.fig_qxz = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Shear XZ component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.qxz_nodeextrap_min mdl.res.qxz_nodeextrap_max]);
                else
                    caxis([mdl.res.qxz_elemextrap_min mdl.res.qxz_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Shear YZ
            if (mdl.res.qyz)
                this.fig_qyz = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Shear YZ component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.qyz_nodeextrap_min mdl.res.qyz_nodeextrap_max]);
                else
                    caxis([mdl.res.qyz_elemextrap_min mdl.res.qyz_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Moment x
            if (mdl.res.mxx)
                this.fig_mxx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Moment XX component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.mxx_nodeextrap_min mdl.res.mxx_nodeextrap_max]);
                else
                    caxis([mdl.res.mxx_elemextrap_min mdl.res.mxx_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Moment Y
            if (mdl.res.myy)
                this.fig_myy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Moment YY component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.myy_nodeextrap_min mdl.res.myy_nodeextrap_max]);
                else
                    caxis([mdl.res.myy_elemextrap_min mdl.res.myy_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Moment XY
            if (mdl.res.mxy)
                this.fig_mxy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Moment XY component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.mxy_nodeextrap_min mdl.res.mxy_nodeextrap_max]);
                else
                    caxis([mdl.res.mxy_elemextrap_min mdl.res.mxy_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Moment 1
            if (mdl.res.m1)
                this.fig_m1 = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum principal moment');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.m1_nodeextrap_min mdl.res.m1_nodeextrap_max]);
                else
                    caxis([mdl.res.m1_elemextrap_min mdl.res.m1_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Moment 2
            if (mdl.res.m2)
                this.fig_m2 = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Minimum principal moment');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.m2_nodeextrap_min mdl.res.m2_nodeextrap_max]);
                else
                    caxis([mdl.res.m2_elemextrap_min mdl.res.m2_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Torsion max
            if (mdl.res.tormax)
                this.fig_tormax = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum torsion moment');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                if (mdl.res.smooth)
                    caxis([mdl.res.tormax_nodeextrap_min mdl.res.tormax_nodeextrap_max]);
                else
                    caxis([mdl.res.tormax_elemextrap_min mdl.res.tormax_elemextrap_max]);
                end
                colorbar;
                hold on;
            end
            
            % Principal moment directions
            if (mdl.res.m1 || mdl.res.m2)
                this.fig_prcmom = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Principal moment directions');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot mesh in active figure.
        function plotMesh(this,mdl)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = this.x_coord(node);
                    YY(j) = this.y_coord(node);
                end
                node1 = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = this.x_coord(node1);
                YY(nen+1) = this.y_coord(node1);
                plot(XX,YY,this.color_mesh);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot mesh ID bumbers for elements, nodes, and Gauss points in
        % active figure.
        function plotMeshLabels(~,mdl)
            if (mdl.res.eid)
                for i = 1:mdl.nel
                    nen = mdl.elems(i).shape.nen;
                    if (nen > 4)
                        nen = nen/2; % ugly workaround for quadratic elements
                    end
                    coords = zeros(nen,2);
                    for j = 1:nen
                        x = mdl.elems(i).shape.nodes(j).coord(1);
                        y = mdl.elems(i).shape.nodes(j).coord(2);
                        coords(j,:) = [x,y];
                    end
                    pol = polyshape;
                    pol.Vertices = coords;
                    [x,y] = centroid(pol);
                    id = mdl.elems(i).id;
                    text(x,y,num2str(id),'color','r');
                    hold on;
                end
            end
            if (mdl.res.nid)
                for i = 1:mdl.nnp
                    id = mdl.nodes(i).id;
                    x = mdl.nodes(i).coord(1);
                    y = mdl.nodes(i).coord(2);
                    scatter(x,y,'b','.');
                    text(x,y,num2str(id),'color','b');
                    hold on;
                end
            end
            if (mdl.res.gid)
                k = 0;
                for i = 1:mdl.nel
                    for j = 1:mdl.res.ngp(i)
                        k = k + 1;
                        x = mdl.res.x_gp(k);
                        y = mdl.res.y_gp(k);
                        scatter(x,y,'g','*');
                        hold on;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plot deformed mesh of 2D inplane analysis model in active figure.
        function plotDeformMeshInplane(this,mdl)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            
            u = mdl.res.U(mdl.ID(1,:));
            v = mdl.res.U(mdl.ID(2,:));
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = this.x_coord(node) + this.deform_fac*u(node);
                    YY(j) = this.y_coord(node) + this.deform_fac*v(node);
                end
                node1 = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = this.x_coord(node1) + this.deform_fac*u(node1);
                YY(nen+1) = this.y_coord(node1) + this.deform_fac*v(node1);
                plot(XX,YY,this.color_deform);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot deformed mesh of 2D outplane analysis model in active figure.
        function plotDeformMeshOutplane(this,mdl)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            w = mdl.res.U(mdl.ID(1,:));
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = this.x_coord(node);
                    YY(j) = this.y_coord(node);
                    ZZ(j) = this.deform_fac*w(node);
                end
                node1 = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = this.x_coord(node1);
                YY(nen+1) = this.y_coord(node1);
                ZZ(nen+1) = this.deform_fac*w(node1);
                plot3(XX,YY,ZZ,this.color_deform);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot element contour response in active figure.
        % Element contour is based on non average element results.
        function plotElemContour(this,mdl,contour_type)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            switch contour_type
                case this.SXX_ELEMEXTRAP
                    contour = mdl.res.sxx_elemextrap;
                case this.SYY_ELEMEXTRAP
                    contour = mdl.res.syy_elemextrap;
                case this.TXY_ELEMEXTRAP
                    contour = mdl.res.txy_elemextrap;
                case this.S1_ELEMEXTRAP
                    contour = mdl.res.s1_elemextrap;
                case this.S2_ELEMEXTRAP
                    contour = mdl.res.s2_elemextrap;
                case this.TMAX_ELEMEXTRAP
                    contour = mdl.res.tmax_elemextrap;
                case this.QXZ_ELEMEXTRAP
                    contour = mdl.res.qxz_elemextrap;
                case this.QYZ_ELEMEXTRAP
                    contour = mdl.res.qyz_elemextrap;
                case this.MXX_ELEMEXTRAP
                    contour = mdl.res.mxx_elemextrap;
                case this.MYY_ELEMEXTRAP
                    contour = mdl.res.myy_elemextrap;
                case this.MXY_ELEMEXTRAP
                    contour = mdl.res.mxy_elemextrap;
                case this.M1_ELEMEXTRAP
                    contour = mdl.res.m1_elemextrap;
                case this.M2_ELEMEXTRAP
                    contour = mdl.res.m2_elemextrap;
                case this.TORMAX_ELEMEXTRAP
                    contour = mdl.res.tormax_elemextrap;
                case this.FXX_ELEMEXTRAP
                    contour = mdl.res.fxx_elemextrap;
                case this.FYY_ELEMEXTRAP
                    contour = mdl.res.fyy_elemextrap;
                case this.FM_ELEMEXTRAP
                    contour = mdl.res.fm_elemextrap;
            end
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node  = mdl.elems(i).shape.ccwNodeIds(j);
                    locid = mdl.elems(i).shape.ccwLocalNodeIds(j);
                    XX(j) = this.x_coord(node);
                    YY(j) = this.y_coord(node);
                    ZZ(j) = contour(locid,i);
                end
                node = mdl.elems(i).shape.ccwNodeIds(1);
                locid = mdl.elems(i).shape.ccwLocalNodeIds(1);
                XX(nen+1) = this.x_coord(node);
                YY(nen+1) = this.y_coord(node);
                ZZ(nen+1) = contour(locid,i);
                patch(XX,YY,ZZ);
                hold on;
            end
            
            % Set callback function to manage nodal data tip showing with cursor
            dcm = datacursormode;
            dcm.UpdateFcn = this.fcn_dataTip;
            dcm.Enable = 'off';
        end
        
        %------------------------------------------------------------------
        % Plot node contour response in active figure.
        % Node contour is based on nodal average adjacent element results.
        function plotNodeContour(this,mdl,contour_type)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
            switch contour_type
                case this.DX_NODE
                    contour = mdl.res.U(mdl.ID(1,:));
                case this.DY_NODE
                    contour = mdl.res.U(mdl.ID(2,:));
                case this.DZ_NODE
                    contour = -mdl.res.U(mdl.ID(1,:));  % THIS IS EXCLUSIVE TO THICK PALTE !!!
                case this.RX_NODE
                    contour = mdl.res.U(mdl.ID(2,:));   % THIS IS EXCLUSIVE TO THICK PALTE !!!
                case this.RY_NODE
                    contour = mdl.res.U(mdl.ID(3,:));   % THIS IS EXCLUSIVE TO THICK PALTE !!!
                case this.SXX_NODEEXTRAP
                    contour = mdl.res.sxx_nodeextrap;
                case this.SYY_NODEEXTRAP
                    contour = mdl.res.syy_nodeextrap;
                case this.TXY_NODEEXTRAP
                    contour = mdl.res.txy_nodeextrap;
                case this.S1_NODEEXTRAP
                    contour = mdl.res.s1_nodeextrap;
                case this.S2_NODEEXTRAP
                    contour = mdl.res.s2_nodeextrap;
                case this.TMAX_NODEEXTRAP
                    contour = mdl.res.tmax_nodeextrap;
                case this.QXZ_NODEEXTRAP
                    contour = mdl.res.qxz_nodeextrap;
                case this.QYZ_NODEEXTRAP
                    contour = mdl.res.qyz_nodeextrap;
                case this.MXX_NODEEXTRAP
                    contour = mdl.res.mxx_nodeextrap;
                case this.MYY_NODEEXTRAP
                    contour = mdl.res.myy_nodeextrap;
                case this.MXY_NODEEXTRAP
                    contour = mdl.res.mxy_nodeextrap;
                case this.M1_NODEEXTRAP
                    contour = mdl.res.m1_nodeextrap;
                case this.M2_NODEEXTRAP
                    contour = mdl.res.m2_nodeextrap;
                case this.TORMAX_NODEEXTRAP
                    contour = mdl.res.tormax_nodeextrap;
                case this.TEMP_NODE
                    contour = mdl.res.U(mdl.ID(1,:));
                case this.FXX_NODEEXTRAP
                    contour = mdl.res.fxx_nodeextrap;
                case this.FYY_NODEEXTRAP
                    contour = mdl.res.fyy_nodeextrap;
                case this.FM_NODEEXTRAP
                    contour = mdl.res.fm_nodeextrap;
            end
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    node  = mdl.elems(i).shape.ccwNodeIds(j);
                    XX(j) = this.x_coord(node);
                    YY(j) = this.y_coord(node);
                    ZZ(j) = contour(node);
                end
                node = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = this.x_coord(node);
                YY(nen+1) = this.y_coord(node);
                ZZ(nen+1) = contour(node);
                patch(XX,YY,ZZ);
                hold on;
            end
            
            % Set callback function to manage nodal data tip showing with cursor
            dcm = datacursormode;
            dcm.UpdateFcn = this.fcn_dataTip;
            dcm.Enable = 'off';
        end
        
        %------------------------------------------------------------------
        % Create figures for plotting results of static analysis.
        function createStaticFigs(this,mdl)
            if (mdl.anm.type == fem.Anm.PLANE_STRESS      || ...
                mdl.anm.type == fem.Anm.PLANE_STRAIN      || ...
                mdl.anm.type == fem.Anm.AXISYM_STRESS     || ...
                mdl.anm.type == fem.Anm.PLANE_CONDUCTION  || ...
                mdl.anm.type == fem.Anm.AXISYM_CONDUCTION || ...
                mdl.anm.type == fem.Anm.CONVECTION_DIFFUSION)
                this.createStaticFigsInplane(mdl);
            elseif (mdl.anm.type == fem.Anm.THICK_PLATE)
                this.createStaticFigsOutplane(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in current active figure.
        function plotDeformMesh(this,mdl)
            if (mdl.anm.type == fem.Anm.PLANE_STRESS || ...
                mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
                mdl.anm.type == fem.Anm.AXISYM_STRESS)
                this.plotDeformMeshInplane(mdl);
            elseif (mdl.anm.type == fem.Anm.THICK_PLATE)
                this.plotDeformMeshOutplane(mdl);
            end
        end
    end
end