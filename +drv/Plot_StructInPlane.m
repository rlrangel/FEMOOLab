%% Plot_StructInPlane Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <plot.html Plot: plot super-class> to deal with
% in-plane plotting of 2D structural models.
%
%% Class definition
%
classdef Plot_StructInPlane < drv.Plot
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Handles to plot figures
        fig_lbl    = [];  % figure for mesh labels plot
        fig_deform = [];  % figure for deformed mesh plot
        fig_dx     = [];  % figure for displacement x plot
        fig_dy     = [];  % figure for displacement y plot
        fig_sxx    = [];  % figure for sigma x plot
        fig_syy    = [];  % figure for sigma y plot
        fig_txy    = [];  % figure for tau xy plot
        fig_s1     = [];  % figure for sigma 1 plot
        fig_s2     = [];  % figure for sigma 2 plot
        fig_tmax   = [];  % figure for tau max. plot
        fig_prcstr = [];  % figure for principal stress directions
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Plot_StructInPlane()
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
            
%             % Deformed configuration
%             if (mdl.res.deform)
%                 % Create nodal displacement response data
%                 udispl     = mdl.res.U(mdl.ID(1,:));
%                 vdispl     = mdl.res.U(mdl.ID(2,:));
%                 udispl_min = min(udispl);
%                 udispl_max = max(udispl);
%                 vdispl_min = min(vdispl);
%                 vdispl_max = max(vdispl);
%                 
%                 % Compute deformed factor based on maximum displacement abs. value
%                 max_displ = max([abs(udispl_min), abs(udispl_max), abs(vdispl_min), abs(vdispl_max)]);
%                 min_plotsize = min(this.plot_xmax - this.plot_xmin, this.plot_ymax - this.plot_ymin);
%                 
%                 if (mdl.res.scl == 0)
%                     this.deform_fac = (min_plotsize / max_displ) * 0.15;
%                 else
%                     this.deform_fac = mdl.res.scl;
%                 end
%                 
%                 % Create figure for mesh and deformed mesh plot
%                 this.fig_deform = figure;
%                 movegui(gca,'center');
%                 set(gca,'DataAspectRatio',[1 1 1]);
%                 title_text = sprintf('Mesh and deformed mesh. Deformed factor: %s',num2str(this.deform_fac));
%                 title(title_text);
%                 shift = max_displ * this.deform_fac;
%                 plot_xmin_shift = this.plot_xmin - shift;
%                 plot_xmax_shift = this.plot_xmax + shift;
%                 plot_ymin_shift = this.plot_ymin - shift;
%                 plot_ymax_shift = this.plot_ymax + shift;
%                 axis([plot_xmin_shift plot_xmax_shift plot_ymin_shift plot_ymax_shift]);
%                 hold on;
%             end
%             
%             % Displacement X
%             if (mdl.res.dx)
%                 this.fig_dx = figure;
%                 movegui(gca,'center');
%                 set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
%                 title('Displacement X');
%                 axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
%                 colorbar;
%                 hold on;
%             end
%             
%             % Displacement Y
%             if (mdl.res.dy)
%                 this.fig_dy = figure;
%                 movegui(gca,'center');
%                 set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
%                 title('Displacement Y');
%                 axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
%                 colorbar;
%                 hold on;
%             end
            
            % Sigma X
            if (mdl.res.sxx)
                this.fig_sxx = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Sigma XX stress component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Sigma Y
            if (mdl.res.syy)
                this.fig_syy = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Sigma YY stress component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Tau XY
            if (mdl.res.txy)
                this.fig_txy = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Tau XY stress component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % sigma 1
            if (mdl.res.s1)
                this.fig_s1 = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum principal stress');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % sigma 2
            if (mdl.res.s2)
                this.fig_s2 = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Minimum principal stress');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Tau max
            if (mdl.res.taumax)
                this.fig_tmax = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Maximum shear stress');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Principal stress directions
            if (mdl.res.s1 || mdl.res.s2)
                this.fig_prcstr = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Principal stress directions');
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
            
%             if (mdl.res.deform)
%                 figure(this.fig_deform);
%                 this.plotMesh2D(mdl);
%                 this.plotDeformMesh(mdl);
%             end
%             
%             if (mdl.res.dx)
%                 figure(this.fig_dx);
%                 contour = mdl.res.U(mdl.ID(1,:));
%                 this.plotNodeContour2D(mdl,contour);
%             end
%             
%             if (mdl.res.dy)
%                 figure(this.fig_dy);
%                 contour = mdl.res.U(mdl.ID(2,:));
%                 this.plotNodeContour2D(mdl,contour);
%             end
            
            if (mdl.res.sxx)
                figure(this.fig_sxx);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.sxx_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.sxx_elemextrap);
                end
            end
            
            if (mdl.res.syy)
                figure(this.fig_syy);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.syy_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.syy_elemextrap);
                end
            end
            
            if (mdl.res.txy)
                figure(this.fig_txy);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.txy_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.txy_elemextrap);
                end
            end
            
            if (mdl.res.s1)
                figure(this.fig_s1);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.s1_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.s1_elemextrap);
                end
            end
            
            if (mdl.res.s2)
                figure(this.fig_s2);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.s2_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.s2_elemextrap);
                end
            end
            
            if (mdl.res.taumax)
                figure(this.fig_tmax);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.tmax_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.tmax_elemextrap);
                end
            end
            
            if (mdl.res.s1 || mdl.res.s2)
                figure(this.fig_prcstr);
                this.plotMesh2D(mdl);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.s1x_gp,mdl.res.s1y_gp,'r');
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.s2x_gp,mdl.res.s2y_gp,'b');
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
%             if (mdl.res.dy)
%                 figure(this.fig_dy);
%             end
%             if (mdl.res.dx)
%                 figure(this.fig_dx);
%             end
%             if (mdl.res.deform)
%                 figure(this.fig_deform);
%             end
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
            u = mdl.res.U(mdl.ID(1,:));
            v = mdl.res.U(mdl.ID(2,:));
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                XX = zeros(1,nen+1);
                YY = zeros(1,nen+1);
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
    end
end