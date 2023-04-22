%% Plot_ThermalPlane Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <plot.html Plot: plot super-class> to deal with
% in-plane plotting of 2D thermal models.
%
%% Class definition
%
classdef Plot_ThermalPlane < drv.Plot
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Handles to plot figures
        fig_lbl      = [];  % figure for mesh labels
        fig_pec      = [];  % figure for element Peclet numbers
        fig_temp     = [];  % figure for temperature field
        fig_fxx      = [];  % figure for flux x plot
        fig_fyy      = [];  % figure for flux x plot
        fig_fm       = [];  % figure for flux module plot
        fig_prcflx   = [];  % figure for flux directions
        fig_crvTemp  = [];  % figure for temperature curve graph
        fig_crvTempt = [];  % figure for temperature rate of change curve graph
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Plot_ThermalPlane()
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
            
            % Peclet number
            if (mdl.res.pec)
                this.fig_pec = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',parula);
                title('Element Peclet number');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Temperature field
            if (mdl.res.temp)
                this.fig_temp = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Temperature field');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Heat flux X
            if (mdl.res.fxx)
                this.fig_fxx = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Heat flux X component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Heat flux Y
            if (mdl.res.fyy)
                this.fig_fyy = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Heat flux Y component');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Heat flux module and directions
            if (mdl.res.fm)
                % Heat flux module
                this.fig_fm = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Heat flux module');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
                
                % Heat flux directions
                this.fig_prcflx = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Heat flux directions');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Create figures (windows) for displaying transient results.
        function createFigsTransient(this,mdl)
            % Mesh labels
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                this.fig_lbl = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1]);
                title('Mesh labels');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                hold on;
            end
            
            % Peclet number
            if (mdl.res.pec)
                this.fig_pec = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',parula);
                title('Element Peclet number');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                colorbar;
                hold on;
            end
            
            % Temperature field
            if (mdl.res.temp)
                this.fig_temp = figure;
                movegui(gca,'center');
                set(gca,'DataAspectRatio',[1 1 1],'Colormap',jet);
                title('Temperature field');
                axis([this.plot_xmin this.plot_xmax this.plot_ymin this.plot_ymax]);
                temp_min = min(min(mdl.res.U));
                temp_max = max(max(mdl.res.U));
                if (temp_min < temp_max)
                    caxis([temp_min temp_max]);
                end
                colorbar;
            end
            
            % Nodal temperature curves
            if (~isempty(mdl.res.curve_temp))
                % Temperature
                this.fig_crvTemp = figure;
                movegui(gca,'center');
                title('Nodal Temperature');
                xlabel('Time');
                ylabel('Temperature');
                grid on;
                hold on;
                
                % Temperature rate of change
                this.fig_crvTempt = figure;
                movegui(gca,'center');
                title('Nodal Temperature Rate of Change');
                xlabel('Time');
                ylabel('Temperature Rate of Change');
                grid on;
                hold on; 
            end
        end
        
        %------------------------------------------------------------------
        % Plot contours of steady state analysis results.
        function plotSteadyStateContours(this,mdl)
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
                this.plotMesh2D(mdl);
                this.plotMeshLabels2D(mdl);
            end
            
            if (mdl.res.pec)
                figure(this.fig_pec);
                if (isempty(mdl.res.maxNen))
                    maxNen = mdl.maxNumElemExtNodes();
                else
                    maxNen = mdl.res.maxNen;
                end
                contour = zeros(maxNen,mdl.nel);
                for i = 1:mdl.nel
                    contour(:,i) = mdl.elems(i).peclet;
                end
                this.plotElemContour2D(mdl,contour);
            end
            
            if (mdl.res.temp)
                figure(this.fig_temp);
                contour = mdl.res.U(mdl.ID(1,:));
                this.plotNodeContour2D(mdl,contour);
            end
            
            if (mdl.res.fxx)
                figure(this.fig_fxx);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.fxx_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.fxx_elemextrap);
                end
            end
            
            if (mdl.res.fyy)
                figure(this.fig_fyy);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.fyy_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.fyy_elemextrap);
                end
            end
            
            if (mdl.res.fm)
                figure(this.fig_fm);
                if (mdl.res.smooth)
                    this.plotNodeContour2D(mdl,mdl.res.fm_nodeextrap);
                else
                    this.plotElemContour2D(mdl,mdl.res.fm_elemextrap);
                end
            end
            
            if (mdl.res.fm)
                figure(this.fig_prcflx);
                quiver(mdl.res.x_gp,mdl.res.y_gp,mdl.res.fmx_gp,mdl.res.fmy_gp,'r');
                this.plotMesh2D(mdl);
            end
            
            % Set order to display (reverse)
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
            if (mdl.res.pec)
                figure(this.fig_pec);
            end
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
            end
        end
        
        %------------------------------------------------------------------
        % Plot contours of transient analysis results.
        function plotTransientContours(this,mdl,anl)
            if (mdl.res.eid || mdl.res.nid || mdl.res.gid)
                figure(this.fig_lbl);
                this.plotMesh2D(mdl);
                this.plotMeshLabels2D(mdl);
            end
            
            if (mdl.res.pec)
                figure(this.fig_pec);
                if (isempty(mdl.res.maxNen))
                    maxNen = mdl.maxNumElemExtNodes();
                else
                    maxNen = mdl.res.maxNen;
                end
                contour = zeros(maxNen,mdl.nel);
                for i = 1:mdl.nel
                    contour(:,i) = mdl.elems(i).peclet;
                end
                this.plotElemContour2D(mdl,contour);
            end
            
            if (mdl.res.temp)
                % Number of frames
                steps = 1:mdl.res.steps;
                steps = steps(rem(steps,mdl.res.output_freq)==0);
                loops = length(steps);
                
                % Create arrays with data to be ploted in all steps
                if (isempty(mdl.res.maxNen))
                    maxNen = mdl.maxNumElemExtNodes();
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
                
                % Plot animation
                for i = 1:loops
                    figure(this.fig_temp);
                    for j = 1:mdl.nel
                        patch(XX(j,:,i),YY(j,:,i),ZZ(j,:,i));
                    end
                    title(['Temperature, t = ',num2str(anl.incr*steps(i))],'FontSize',12)
                    pause(eps);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plot graph curves of transient analysis results.
        function plotTransientCurves(this,mdl)
            % Nodal temperature
            if (~isempty(mdl.res.curve_temp))
                ncurves = length(mdl.res.curve_temp);
                for i = 1:ncurves
                    dof = mdl.ID(1,mdl.res.curve_temp(i));
                    x   = mdl.res.times;
                    y   = mdl.res.U(dof,:);
                    yt  = mdl.res.Ut(dof,:);
                    leg = sprintf('Node %d',mdl.res.curve_temp(i));
                    
                    % Plot temperature values
                    figure(this.fig_crvTemp);
                    plot(x,y,'DisplayName',leg);
                    
                    % Plot temperature rate of cahnge
                    figure(this.fig_crvTempt);
                    plot(x,yt,'DisplayName',leg);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in active figure.
        function plotDeformMesh(~,~)
            return;
        end
    end
end