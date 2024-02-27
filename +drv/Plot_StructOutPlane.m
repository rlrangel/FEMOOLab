%% Plot_StructOutPlane Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <plot.html Plot: plot super-class> to deal with
% out-of-plane plotting of 2D structural models.
%
%% Class definition
%
classdef Plot_StructOutPlane < drv.Plot
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Handles to plot figures
        fig_lbl    = [];  % figure for mesh labels plot
        fig_deform = [];  % figure for deformed mesh plot
        fig_dz     = [];  % figure for displacement z plot
        fig_rx     = [];  % figure for rotation x plot
        fig_ry     = [];  % figure for rotation y plot
        fig_qxz    = [];  % figure for shear xz plot
        fig_qyz    = [];  % figure for shear yz plot
        fig_mxx    = [];  % figure for moment x plot
        fig_myy    = [];  % figure for moment y plot
        fig_mxy    = [];  % figure for moment xy plot
        fig_m1     = [];  % figure for moment 1 plot
        fig_m2     = [];  % figure for moment 2 plot
        fig_tormax = [];  % figure for torsion max plot
        fig_prcmom = [];  % figure for principal moment directions
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Plot_StructOutPlane()
            this = this@drv.Plot();
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Plot
    methods
        %------------------------------------------------------------------
        % Create figures (windows) for displaying steady state results.
        function createFigsSteadyState(this,mdl)
            % Mesh labels
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                this.fig_lbl = figure;
                movegui(gca,'center');
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
                movegui(gca,'center');
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
                this.fig_dz = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Displacement Z');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Rotation X
            if (mdl.res.rx)
                this.fig_rx = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Rotation X');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Rotation Y
            if (mdl.res.ry)
                this.fig_ry = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Rotation Y');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Shear XZ
            if (mdl.res.qxz)
                this.fig_qxz = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Shear XZ component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Shear YZ
            if (mdl.res.qyz)
                this.fig_qyz = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Shear YZ component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Moment x
            if (mdl.res.mxx)
                this.fig_mxx = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Moment XX component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Moment Y
            if (mdl.res.myy)
                this.fig_myy = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Moment YY component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Moment XY
            if (mdl.res.mxy)
                this.fig_mxy = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Moment XY component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Moment 1
            if (mdl.res.m1)
                this.fig_m1 = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum principal moment');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Moment 2
            if (mdl.res.m2)
                this.fig_m2 = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Minimum principal moment');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Torsion max
            if (mdl.res.tormax)
                this.fig_tormax = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum torsion moment');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Principal moment directions
            if (mdl.res.m1 || mdl.res.m2)
                this.fig_prcmom = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Principal moment directions');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Create figures (windows) for displaying transient results.
        function createFigsTransient(~,~)
            % NOT IMPLEMENTED !!!
            return;
        end
        
        %------------------------------------------------------------------
        % Plot contours of steady state analysis results.
        function plotSteadyStateContours(this,mdl)
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
                this.plotMesh2D(mdl);
                this.plotMeshLabels2D(mdl);
            end
            
            if (mdl.res.deform)
                figure(this.fig_deform);
                this.plotMesh2D(mdl);
                this.plotDeformMesh(mdl);
            end
            
            if (mdl.res.dz)
                figure(this.fig_dz);
                contour = -mdl.res.U(mdl.ID(1,:));
                this.plotNodeContour2D(mdl,contour);
            end
            
            if (mdl.res.rx)
                figure(this.fig_rx);
                contour = mdl.res.U(mdl.ID(2,:));
                this.plotNodeContour2D(mdl,contour);
            end
            
            if (mdl.res.ry)
                figure(this.fig_ry);
                contour = mdl.res.U(mdl.ID(3,:));
                this.plotNodeContour2D(mdl,contour);
            end
            
            if (mdl.res.qxz)
                figure(this.fig_qxz);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.qxz_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.qxz_elemextrap);
                end
            end
            
            if (mdl.res.qyz)
                figure(this.fig_qyz);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.qyz_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.qyz_elemextrap);
                end
            end
            
            if (mdl.res.mxx)
                figure(this.fig_mxx);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.mxx_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.mxx_elemextrap);
                end
            end
            
            if (mdl.res.myy)
                figure(this.fig_myy);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.myy_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.myy_elemextrap);
                end
            end
            
            if (mdl.res.mxy)
                figure(this.fig_mxy);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.mxy_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.mxy_elemextrap);
                end
            end
            
            if (mdl.res.m1)
                figure(this.fig_m1);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.m1_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.m1_elemextrap);
                end
            end
            
            if (mdl.res.m2)
                figure(this.fig_m2);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.m2_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.m2_elemextrap);
                end
            end
            
            if (mdl.res.tormax)
                figure(this.fig_tormax);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.tormax_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.tormax_elemextrap);
                end
            end
            
            if (mdl.res.m1 || mdl.res.m2)
                figure(this.fig_prcmom);
                this.plotMesh2D(mdl);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.m1x_gp,mdl.res.m1y_gp,'r');
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.m2x_gp,mdl.res.m2y_gp,'b');
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
            if (mdl.res.ry)
                figure(this.fig_ry);
            end
            if (mdl.res.rx)
                figure(this.fig_rx);
            end
            if (mdl.res.dz)
                figure(this.fig_dz);
            end
            if (mdl.res.deform)
                figure(this.fig_deform);
            end
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
            end
        end
        
        %------------------------------------------------------------------
        % Plot contours of transient analysis results.
        function plotTransientContours(~,~,~)
            % NOT IMPLEMENTED !!!
            return;
        end
        
        %------------------------------------------------------------------
        % Plot graph curves of transient analysis results.
        function plotTransientCurves(~,~)
            % NOT IMPLEMENTED !!!
            return;
        end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in active figure.
        function plotDeformMesh(this,mdl)
            w = mdl.res.U(mdl.ID(1,:));
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                XX = zeros(1,nen+1);
                YY = zeros(1,nen+1);
                ZZ = zeros(1,nen+1);
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
    end
end