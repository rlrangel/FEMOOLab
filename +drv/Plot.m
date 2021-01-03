%% Plot Class
%
%% Description
%
% This class defines a plotting object in the FEMOOLab program.
% A plotting object is responsible for displaying the analysis results in
% Matlab figure windows.
%
classdef Plot < handle
    %% Constant values for contour types
    properties (Constant = true, Access = public)
        % Structural analysis
        DX_NODE         = int32(1);
        DY_NODE         = int32(2);
        SXX_ELEMEXTRAP  = int32(4);
        SYY_ELEMEXTRAP  = int32(5);
        TXY_ELEMEXTRAP  = int32(6);
        S1_ELEMEXTRAP   = int32(7);
        S2_ELEMEXTRAP   = int32(8);
        TMAX_ELEMEXTRAP = int32(9);
        SXX_NODEEXTRAP  = int32(10);
        SYY_NODEEXTRAP  = int32(11);
        TXY_NODEEXTRAP  = int32(12);
        S1_NODEEXTRAP   = int32(13);
        S2_NODEEXTRAP   = int32(14);
        TMAX_NODEEXTRAP = int32(15);
        
        % Thermal analysis
        TEMP_NODE       = int32(16);
        FXX_ELEMEXTRAP  = int32(17);
        FYY_ELEMEXTRAP  = int32(18);
        FM_ELEMEXTRAP   = int32(19);
        FXX_NODEEXTRAP  = int32(20);
        FYY_NODEEXTRAP  = int32(21);
        FM_NODEEXTRAP   = int32(22);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Handles to plot figures
        fig_lbl    = [];  % figure for mesh labels
        fig_deform = [];  % figure for mesh and deformed mesh plot
        fig_dx     = [];  % figure for displacement x plot
        fig_dy     = [];  % figure for displacement y plot
        fig_sxx    = [];  % figure for sigma x plot
        fig_syy    = [];  % figure for sigma y plot
        fig_txy    = [];  % figure for tau xy plot
        fig_s1     = [];  % figure for sigma 1 plot
        fig_s2     = [];  % figure for sigma 2 plot
        fig_tmax   = [];  % figure for tau max. plot
        fig_prcstr = [];  % figure for principal stress directions
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
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Plot(sim)
            if (nargin > 0)
                % Setup bounding box for displaying results
                this.setupBoundingBox(sim.mdl);
                
                if (sim.anl.type == fem.Anl.LINEAR_STATIC)
                    this.plotStatic(sim.mdl);
                elseif (sim.anl.type == fem.Anl.LINEAR_TRANSIENT)
                    
                end
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Plot static analysis results.
        function plotStatic(this,mdl)
            if (mdl.anm.type == fem.Anm.PLANE_STRESS     || ...
                mdl.anm.type == fem.Anm.PLANE_STRAIN     || ...
                mdl.anm.type == fem.Anm.AXISYM_STRESS    || ...
                mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                mdl.anm.type == fem.Anm.AXISYM_CONDUCTION)
                this.plotStaticInplane(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Create figures for plotting results of static analysis.
        function deform_fac = createStaticFigs(this,mdl)
            if (mdl.anm.type == fem.Anm.PLANE_STRESS     || ...
                mdl.anm.type == fem.Anm.PLANE_STRAIN     || ...
                mdl.anm.type == fem.Anm.AXISYM_STRESS    || ...
                mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                mdl.anm.type == fem.Anm.AXISYM_CONDUCTION)
                deform_fac = this.createStaticFigsInplane(mdl);
            end
        end
        
        %------------------------------------------------------------------
        % Plot mesh in current active figure.
        function plotMesh(this,mdl,color)
            if (mdl.anm.type == fem.Anm.PLANE_STRESS     || ...
                mdl.anm.type == fem.Anm.PLANE_STRAIN     || ...
                mdl.anm.type == fem.Anm.AXISYM_STRESS    || ...
                mdl.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                mdl.anm.type == fem.Anm.AXISYM_CONDUCTION)
                this.plotMeshInplane(mdl,color);
            end
       end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in current active figure.
        function plotDeformMesh(this,mdl,fct,color)
            if (mdl.anm.type == fem.Anm.PLANE_STRESS || ...
                mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
                mdl.anm.type == fem.Anm.AXISYM_STRESS)
                this.plotDeformMeshInplane(mdl,fct,color);
            end
        end
    end
    
    %% Private methods
    methods (Access = private)
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
        % Plot 2D inplane static analysis results.
        function plotStaticInplane(this,mdl)
            % Create figures (windows) for displaying results
            deform_fac = this.createStaticFigs(mdl);
            
            % Display results
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
                this.plotMesh(mdl,'k');
                this.plotMeshLabels(mdl);
            end
            
            if (mdl.res.deform)
                figure(this.fig_deform);
                this.plotMesh(mdl,'k');
                this.plotDeformMesh(mdl,deform_fac,'b');
            end
            
            if (mdl.res.dx)
                figure(this.fig_dx);
                this.plotNodeContourInplane(mdl,this.DX_NODE);
            end
            
            if (mdl.res.dy)
                figure(this.fig_dy);
                this.plotNodeContourInplane(mdl,this.DY_NODE);
            end
            
            if (mdl.res.sxx)
                figure(this.fig_sxx);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.SXX_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.SXX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.syy)
                figure(this.fig_syy);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.SYY_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.SYY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.txy)
                figure(this.fig_txy);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.TXY_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.TXY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.s1)
                figure(this.fig_s1);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.S1_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.S1_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.s2)
                figure(this.fig_s2);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.S2_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.S2_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.taumax)
                figure(this.fig_tmax);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.TMAX_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.TMAX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.s1 || mdl.res.s2)
                figure(this.fig_prcstr);
                this.plotMesh(mdl,'k');
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.s1x_gp,mdl.res.s1y_gp,'r');
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.s2x_gp,mdl.res.s2y_gp,'b');
            end
            
            if (mdl.res.temp)
                figure(this.fig_temp);
                this.plotNodeContourInplane(mdl,this.TEMP_NODE);
            end
            
            if (mdl.res.fxx)
                figure(this.fig_fxx);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.FXX_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.FXX_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.fyy)
                figure(this.fig_fyy);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.FYY_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.FYY_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.fm)
                figure(this.fig_fm);
                if (mdl.res.smooth)
                    this.plotNodeContourInplane(mdl,this.FM_NODEEXTRAP);
                else
                    this.plotElemContourInplane(mdl,this.FM_ELEMEXTRAP);
                end
            end
            
            if (mdl.res.fm)
                figure(this.fig_prcflx);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.fmx_gp,mdl.res.fmy_gp,'r');
                this.plotMesh(mdl,'k');
            end
            
            % Set order to display (reverse)
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
        % Create figures for plotting 2D inplane static analysis results.
        function deform_fac = createStaticFigsInplane(this,mdl)
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
                    deform_fac = (min_plotsize / max_displ) * 0.15;
                else
                    deform_fac = mdl.res.scl;
                end
                
                % Create figure for mesh and deformed mesh plot and get its handle
                this.fig_deform = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title_text = sprintf('Mesh and deformed mesh. Deformed factor: %s',num2str(deform_fac));
                title(title_text);
                shift = max_displ * deform_fac;
                plot_xmin_shift = this.plot_xmin - shift;
                plot_xmax_shift = this.plot_xmax + shift;
                plot_ymin_shift = this.plot_ymin - shift;
                plot_ymax_shift = this.plot_ymax + shift;
                axis([plot_xmin_shift plot_xmax_shift plot_ymin_shift plot_ymax_shift]);
                hold on;
            else
                deform_fac = 0;
            end
            
            % Displacement X
            if (mdl.res.dx)
                % Create nodal displacement response data
                udispl = mdl.res.U(mdl.ID(1,:));
                udispl_min = min(udispl);
                udispl_max = max(udispl);
                
                % Create figure for temperature plot and get its handle.
                this.fig_dx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Displacement X');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                caxis([udispl_min udispl_max]);
                colorbar;
                hold on;
            end
            
            % Displacement X
            if (mdl.res.dy)
                % Create nodal displacement response data
                vdispl = mdl.res.U(mdl.ID(2,:));
                vdispl_min = min(vdispl);
                vdispl_max = max(vdispl);
                
                % Create figure for temperature plot and get its handle.
                this.fig_dy = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Displacement Y');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                caxis([vdispl_min vdispl_max]);
                colorbar;
                hold on;
            end
            
            % Sigma X
            if (mdl.res.sxx)
                this.fig_sxx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
                % Create nodal temperature response data
                temperature = mdl.res.U(mdl.ID(1,:));
                temp_min = min(temperature);
                temp_max = max(temperature);
                
                % Create figure for temperature plot and get its handle.
                this.fig_temp = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Temperature field');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                caxis([temp_min temp_max]);
                colorbar;
                hold on;
            end
            
            % Heat flux X
            if (mdl.res.fxx)
                this.fig_fxx = figure;
                movegui(gca,'center')
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
                set(gca,'DataAspectRatio',[1 1 1]);
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
        % Plot mesh of 2D inplane analysis model in active figure.
        function plotMeshInplane(this,mdl,color)
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
                plot(XX,YY,color);
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
        function plotDeformMeshInplane(this,mdl,fct,color)
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
                    XX(j) = this.x_coord(node) + fct*u(node);
                    YY(j) = this.y_coord(node) + fct*v(node);
                end
                node1 = mdl.elems(i).shape.ccwNodeIds(1);
                XX(nen+1) = this.x_coord(node1) + fct*u(node1);
                YY(nen+1) = this.y_coord(node1) + fct*v(node1);
                plot(XX,YY,color);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot element contour response of 2D inplane analysis model in
        % active figure.
        % Element contour is based on non average element results.
        function plotElemContourInplane(this,mdl,contour_type)
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
        end
        
        %------------------------------------------------------------------
        % Plot node contour response of 2D inplane analysis model in 
        % active figure.
        % Node contour is based on nodal average adjacent element results.
        function plotNodeContourInplane(this,mdl,contour_type)
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
                    contour = mdl.res.U(mdl.ID(1,1:mdl.nnp));
                case this.DY_NODE
                    contour = mdl.res.U(mdl.ID(2,1:mdl.nnp));
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
                case this.TEMP_NODE
                    contour = mdl.res.U(mdl.ID(1,1:mdl.nnp));
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
        end
    end
end