%% Result Class
%
% This class defines a result object in the StAnOOP program.
% A Result object is responsible for storing the analysis results and
% provide them in text or graphical formats.
%
classdef Result < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        D = [];             % global displacement vector
        
        % Gauss point results
        x_gp;               % vector of gauss point x coordinates
        y_gp;               % vector of gauss point y coordinates
        sx_gp               % sigma x gauss points stress array
        sx_gp_min           % sigma x gauss points min. stress
        sx_gp_max           % sigma x gauss points max. stress
        sy_gp               % sigma y gauss points stress array
        sy_gp_min           % sigma y gauss points min. stress
        sy_gp_max           % sigma y gauss points max. stress
        txy_gp              % tau xy gauss points stress array
        txy_gp_min          % tau xy gauss points min. stress
        txy_gp_max          % tau xy gauss points max. stress
        s1_gp               % sigma 1 gauss points stress array
        s1_gp_min           % sigma 1 gauss points min. stress
        s1_gp_max           % sigma 1 gauss points max. stress
        s2_gp               % sigma 2 gauss points stress array
        s2_gp_min           % sigma 2 gauss points min. stress
        s2_gp_max           % sigma 2 gauss points max. stress
        tmax_gp             % tau max. gauss points stress array
        tmax_gp_min         % tau max. gauss points min. stress
        tmax_gp_max         % tau max. gauss points max. stress
        s1x_gp              % vector of s1 vector gauss x components
        s1y_gp              % vector of s1 vector gauss y components
        s2x_gp              % vector of s2 vector gauss x components
        s2y_gp              % vector of s2 vector gauss y components
        
        % Element node results
        sx_elemextrap       % sigma x element node extrap. stress array
        sx_elemextrap_min   % sigma x element node extrap. min. stress
        sx_elemextrap_max   % sigma x element node extrap. max. stress
        sy_elemextrap       % sigma y element node extrap. stress array
        sy_elemextrap_min   % sigma y element node extrap. min. stress
        sy_elemextrap_max   % sigma y element node extrap. max. stress
        txy_elemextrap      % tau xy element node extrap. stress array
        txy_elemextrap_min  % tau xy element node extrap. min. stress
        txy_elemextrap_max  % tau xy element node extrap. max. stress
        s1_elemextrap       % sigma 1 element node extrap. stress array
        s1_elemextrap_min   % sigma 1 element node extrap. min. stress
        s1_elemextrap_max   % sigma 1 element node extrap. max. stress
        s2_elemextrap       % sigma 2 element node extrap. stress array
        s2_elemextrap_min   % sigma 2 element node extrap. min. stress
        s2_elemextrap_max   % sigma 2 element node extrap. max. stress
        tmax_elemextrap     % tau max. element node extrap. stress array
        tmax_elemextrap_min % tau max. element node extrap. min stress
        tmax_elemextrap_max % tau max. element node extrap. max stress
        
        % Global node results
        sx_nodeextrap       % sigma x extrap. node smoothed stress array
        sx_nodeextrap_min   % sigma x extrap. node smoothed min. stress
        sx_nodeextrap_max   % sigma x extrap. node smoothed max. stress
        sy_nodeextrap       % sigma y extrap. node smoothed stress array
        sy_nodeextrap_min   % sigma y extrap. node smoothed min. stress
        sy_nodeextrap_max   % sigma y extrap. node smoothed max. stress
        txy_nodeextrap      % tau xy extrap. node smoothed stress array
        txy_nodeextrap_min  % tau xy extrap. node smoothed min. stress
        txy_nodeextrap_max  % tau xy extrap. node smoothed max. stress
        s1_nodeextrap       % sigma 1 extrap. node smoothed stress array
        s1_nodeextrap_min   % sigma 1 extrap. node smoothed min. stress
        s1_nodeextrap_max   % sigma 1 extrap. node smoothed max. stress
        s2_nodeextrap       % sigma 2 extrap. node smoothed stress array
        s2_nodeextrap_min   % sigma 2 extrap. node smoothed min. stress
        s2_nodeextrap_max   % sigma 2 extrap. node smoothed max. stress
        tmax_nodeextrap     % tau max. extrap. node smoothed stress array
        tmax_nodeextrap_min % tau max. extrap. node smoothed min. stress
        tmax_nodeextrap_max % tau max. extrap. node smoothed max. stress
        
        % Handles to plot figures
        fig_deform =[];     % figure for mesh and deformed mesh plot
        fig_strbar = [];    % figure for stress bar response plots
        fig_sx     = [];    % figure for sigma x plot
        fig_sy     = [];    % figure for sigma y plot
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
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Clear numerical garbage from analysis. (TEMPORARY SIMPLIFICATION !!)
        function clearSmallValues(res,mdl)
            if (abs(res.sx_gp_min) < 0.00001 && abs(res.sx_gp_max) < 0.00001)
                res.sx_gp_min = 0.0;
                res.sx_gp_max = 0.0;
                res.sx_gp = zeros(mdl.elems(1).n_gaussstress_pts,mdl.nel);
            end
            if (abs(res.sy_gp_min) < 0.00001 && abs(res.sy_gp_max) < 0.00001)
                res.sy_gp_min = 0.0;
                res.sy_gp_max = 0.0;
                res.sy_gp = zeros(mdl.elems(1).n_gaussstress_pts,mdl.nel);
            end
            if (abs(res.txy_gp_min) < 0.00001 && abs(res.txy_gp_max) < 0.00001)
                res.txy_gp_min = 0.0;
                res.txy_gp_max = 0.0;
                res.txy_gp = zeros(mdl.elems(1).n_gaussstress_pts,mdl.nel);
            end
            if (abs(res.s1_gp_min) < 0.00001 && abs(res.s1_gp_max) < 0.00001)
                res.s1_gp_min = 0.0;
                res.s1_gp_max = 0.0;
                res.s1_gp = zeros(mdl.elems(1).n_gaussstress_pts,mdl.nel);
            end
            if (abs(res.s2_gp_min) < 0.00001 && abs(res.s2_gp_max) < 0.00001)
                res.s2_gp_min = 0.0;
                res.s2_gp_max = 0.0;
                res.s2_gp = zeros(mdl.elems(1).n_gaussstress_pts,mdl.nel);
            end
            if (abs(res.tmax_gp_min) < 0.00001 && abs(res.tmax_gp_max) < 0.00001)
                res.tmax_gp_min = 0.0;
                res.tmax_gp_max = 0.0;
                res.tmax_gp = zeros(mdl.elems(1).n_gaussstress_pts,mdl.nel);
            end
            if (abs(res.sx_elemextrap_min) < 0.00001 && abs(res.sx_elemextrap_max) < 0.00001)
                res.sx_elemextrap_min = 0.0;
                res.sx_elemextrap_max = 0.0;
                res.sx_elemextrap = zeros(mdl.elems(1).nen,mdl.nel);
            end
            if (abs(res.sy_elemextrap_min) < 0.00001 && abs(res.sy_elemextrap_max) < 0.00001)
                res.sy_elemextrap_min = 0.0;
                res.sy_elemextrap_max = 0.0;
                res.sy_elemextrap = zeros(mdl.elems(1).nen,mdl.nel);
            end
            if (abs(res.txy_elemextrap_min) < 0.00001 && abs(res.txy_elemextrap_max) < 0.00001)
                res.txy_elemextrap_min = 0.0;
                res.txy_elemextrap_max = 0.0;
                res.txy_elemextrap = zeros(mdl.elems(1).nen,mdl.nel);
            end
            if (abs(res.s1_elemextrap_min) < 0.00001 && abs(res.s1_elemextrap_max) < 0.00001)
                res.s1_elemextrap_min = 0.0;
                res.s1_elemextrap_max = 0.0;
                res.s1_elemextrap = zeros(mdl.elems(1).nen,mdl.nel);
            end
            if (abs(res.s2_elemextrap_min) < 0.00001 && abs(res.s2_elemextrap_max) < 0.00001)
                res.s2_elemextrap_min = 0.0;
                res.s2_elemextrap_max = 0.0;
                res.s2_elemextrap = zeros(mdl.elems(1).nen,mdl.nel);
            end
            if (abs(res.tmax_elemextrap_min) < 0.00001 && abs(res.tmax_elemextrap_max) < 0.00001)
                res.tmax_elemextrap_min = 0.0;
                res.tmax_elemextrap_max = 0.0;
                res.tmax_elemextrap = zeros(mdl.elems(1).nen,mdl.nel);
            end
        end
        
        %------------------------------------------------------------------
        % Plot analysis results.
        function plot(res,mdl)
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
            figure(res.fig_sx);
            res.plotNodeContour(mdl,x,y,1);
            res.plotMesh(mdl,x,y,'k');
            
            figure(res.fig_sy);
            res.plotNodeContour(mdl,x,y,2);
            res.plotMesh(mdl,x,y,'k');
            
            figure(res.fig_txy);
            res.plotNodeContour(mdl,x,y,3);
            res.plotMesh(mdl,x,y,'k');
            
            figure(res.fig_s1);
            res.plotNodeContour(mdl,x,y,4);
            res.plotMesh(mdl,x,y,'k');
            
            figure(res.fig_s2);
            res.plotNodeContour(mdl,x,y,5);
            res.plotMesh(mdl,x,y,'k');
            
            figure(res.fig_tmax);
            res.plotNodeContour(mdl,x,y,6);
            res.plotMesh(mdl,x,y,'k');
        end
        
        %------------------------------------------------------------------
        % Create figures for post-processing results.
        % Eight windows are created and positioned on the screen.
        % Each window plots a type of post-process result.
        function deform_fac = createFigs(res,mdl,plot_xmin,plot_xmax,plot_ymin,plot_ymax)
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
             res.fig_sx = figure;
             fig_sx_pos = get( res.fig_sx, 'Position' );
             fig_sx_pos(1) = 0;
             fig_sx_pos(2) = (screen_sizes(4) - fig_sx_pos(4))/2;
             set( res.fig_sx, 'Position', fig_sx_pos );
             title( 'Sigma X stress component' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             caxis([res.sx_nodeextrap_min res.sx_nodeextrap_max]);
             colorbar;
             hold on;
            
             % Create figure for sigma y plot and get its handle.
             % Locate figure at the second level right side of screen.
             res.fig_sy = figure;
             fig_sy_pos = get( res.fig_sy, 'Position' );
             fig_sy_pos(1) = screen_sizes(3) - fig_sy_pos(3);
             fig_sy_pos(2) = (screen_sizes(4) - fig_sy_pos(4))/2;
             set( res.fig_sy, 'Position', fig_sy_pos );
             title( 'Sigma Y stress component' );
             set( gca,'DataAspectRatio',[1 1 1] );
             axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
             caxis([res.sy_nodeextrap_min res.sy_nodeextrap_max]);
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
             caxis([res.txy_nodeextrap_min res.txy_nodeextrap_max]);
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
             caxis([res.s1_nodeextrap_min res.s1_nodeextrap_max]);
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
             caxis([res.s2_nodeextrap_min res.s2_nodeextrap_max]);
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
             caxis([res.tmax_nodeextrap_min res.tmax_nodeextrap_max]);
             colorbar;
             hold on;
        end
        
        %------------------------------------------------------------------
        % Plot mesh in current active figure.
        function plotMesh(~,mdl,x,y,color)
            nen = mdl.elems(1).nen;
            
            XX = zeros(1,nen+1);
            YY = zeros(1,nen+1);
            
            % Display mesh
            for i = 1:mdl.nel
                for j = 1:nen
                    node = mdl.elems(i).nodes(j).id;
                    XX(j) = x(node);
                    YY(j) = y(node);
                end
                node1 = mdl.elems(i).nodes(1).id;
                XX(nen+1) = x(node1);
                YY(nen+1) = y(node1);
                plot(XX,YY, color);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot deformed mesh in current active figure.
        function plotDeformMesh(res,mdl,x,y,fct,color)
            nen = mdl.elems(1).nen;
            
            XX = zeros(1,nen+1);
            YY = zeros(1,nen+1);
            
            u = res.D(mdl.ID(1,:));
            v = res.D(mdl.ID(2,:));
            
            % Display deformed mesh
            for i = 1:mdl.nel
                for j = 1:nen
                    node = mdl.elems(i).nodes(j).id;
                    XX(j) = x(node) + fct*u(node);
                    YY(j) = y(node) + fct*v(node);
                end
                node1 = mdl.elems(i).nodes(1).id;
                XX(nen+1) = x(node1) + fct*u(node1);
                YY(nen+1) = y(node1) + fct*v(node1);
                plot(XX,YY,color);
                hold on;
            end
        end
        
        %------------------------------------------------------------------
        % Plot current active node contour response in current active figure.
        function plotNodeContour(res,mdl,x,y,response_type)
            nen = mdl.elems(1).nen;
            
            XX = zeros(1,nen+1);
            YY = zeros(1,nen+1);
            ZZ = zeros(1,nen+1);
            
            % Display current active node contour response
            switch response_type
                case 1
                    contour = res.sx_nodeextrap;
                case 2
                    contour = res.sy_nodeextrap;
                case 3
                    contour = res.txy_nodeextrap;
                case 4
                    contour = res.s1_nodeextrap;
                case 5
                    contour = res.s2_nodeextrap;
                case 6
                    contour = res.tmax_nodeextrap;
            end
            
            for i = 1:mdl.nel
                for j = 1:nen
                    node = mdl.elems(i).nodes(j).id;
                    XX(j) = x(node);
                    YY(j) = y(node);
                    ZZ(j) = contour(node);
                end
                node1 = mdl.elems(i).nodes(1).id;
                XX(nen+1) = x(node1);
                YY(nen+1) = y(node1);
                ZZ(nen+1) = contour(node1);
                patch(XX,YY,ZZ);
                hold on;
            end
        end
    end
end