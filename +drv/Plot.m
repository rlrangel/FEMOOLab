%% Plot Class
%
%% Description
%
% This is an abstract super-class that generically specifies a plotting 
% object object in the FEMOOLab program.
% A plotting object is responsible for displaying the analysis results in
% Matlab figure windows.
% Essentially, this super-class declares abstract methods that are the 
% functions that should be implemented in a derived sub-class that deals
% with specific types of models.
%
%% Subclasses
%
% * <plot_structinplane.html Plot_StructInPlane: structural in-plane plotting subclass>
% * <plot_structoutplane.html Plot_StructOutPlane: structural out-of-plane plotting subclass>
% * <plot_thermalplane.html Plot_ThermalPlane: thermal in-plane plotting subclass>
%
%% Class definition
%
classdef Plot < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Bounding box for plotting figures
        x_coord    double = double.empty;
        y_coord    double = double.empty;
        z_coord    double = double.empty;
        plot_xmin  double = double.empty;
        plot_xmax  double = double.empty;
        plot_ymin  double = double.empty;
        plot_ymax  double = double.empty;
        plot_zmin  double = double.empty;
        plot_zmax  double = double.empty;
        
        % Visual properties
        deform_fac   double = double.empty; % Scale factor for deformed mesh
        color_mesh   = 'k';                 % mesh color
        color_deform = 'b';                 % deformed mesh color
        
        % handle to function to manage data tip showing
        fcn_dataTip  = @drv.cb_dataTipCursor;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Plot()
            return;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Create figures (windows) for displaying steady state results.
        createFigsSteadyState(this,mdl);
        
        %------------------------------------------------------------------
        % Create figures (windows) for displaying transient results.
        createFigsTransient(this,mdl);
        
        %------------------------------------------------------------------
        % Plot contours of steady state analysis results.
        plotSteadyStateContours(this,mdl);
        
        %------------------------------------------------------------------
        % Plot contours of transient analysis results.
        plotTransientContours(this,mdl,anl);
        
        %------------------------------------------------------------------
        % Plot graph curves of transient analysis results.
        plotTransientCurves(this,mdl);
        
        %------------------------------------------------------------------
        % Plot deformed mesh in active figure.
        plotDeformMesh(this,mdl);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Execute plotting of results.
        function execute(this,sim)
            % Set bounding box
            this.setupBoundingBox(sim.mdl);
            
            % Plot results
            if (sim.anl.type == fem.Anl.LINEAR_STEADYSTATE)
                this.createFigsSteadyState(sim.mdl);
                this.plotSteadyStateContours(sim.mdl);
            elseif (sim.anl.type == fem.Anl.LINEAR_TRANSIENT)
                this.createFigsTransient(sim.mdl);
                this.plotTransientCurves(sim.mdl);
                this.plotTransientContours(sim.mdl,sim.anl);
            end
        end
        
        %------------------------------------------------------------------
        % Setup bounding box for plotting figures.
        function setupBoundingBox(this,mdl)
            % Assemble vectors of nodal coordinates
            x = zeros(mdl.nnp,1);
            y = zeros(mdl.nnp,1);
            z = zeros(mdl.nnp,1);
            for i = 1:mdl.nnp
                x(i) = mdl.nodes(i).coord(1);
                y(i) = mdl.nodes(i).coord(2);
                z(i) = mdl.nodes(i).coord(3);
            end
            this.x_coord = x;
            this.y_coord = y;
            this.z_coord = z;
            
            % Setup bounding box for displaying results
            min_x = min(x);
            max_x = max(x);
            min_y = min(y);
            max_y = max(y);
            min_z = min(z);
            max_z = max(z);
            
            size_x = max_x-min_x;
            size_y = max_y-min_y;
            size_z = max_z-min_z;
            
            cx = min_x + size_x * 0.5;
            cy = min_y + size_y * 0.5;
            cz = min_z + size_z * 0.5;
            
            this.plot_xmin = cx - size_x * 0.55;
            this.plot_xmax = cx + size_x * 0.55;
            this.plot_ymin = cy - size_y * 0.55;
            this.plot_ymax = cy + size_y * 0.55;
            this.plot_zmin = cz - size_z * 0.55;
            this.plot_zmax = cz + size_z * 0.55;
        end
        
        %------------------------------------------------------------------
        % Plot 2D mesh in active figure.
        function plotMesh2D(this,mdl)
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
        function plotMeshLabels2D(~,mdl)
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
        % Plot element contour response in active figure.
        % Element contour is based on non average element results.
        function plotElemContour2D(this,mdl,contour)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
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
        function plotNodeContour2D(this,mdl,contour)
            if (isempty(mdl.res.maxNen))
                maxNen = mdl.maxNumElemNodes();
            else
                maxNen = mdl.res.maxNen;
            end
            
            XX = zeros(1,maxNen+1);
            YY = zeros(1,maxNen+1);
            ZZ = zeros(1,maxNen+1);
            
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
    end
end