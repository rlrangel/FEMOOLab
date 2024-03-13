function posCreateFigs
% Creates figures for post-processing results.
% Eight windows are created and positioned on the screen.
% Each window plots a type of post-process result.

 include_gblrefs;

 % Get current screen sizes.
 screen_sizes = get(0,'ScreenSize');

 % Compute deformed factor based on maximum displacement abs. value
 max_displ = max([abs(udispl_min) abs(udispl_max) ...
                  abs(vdispl_min) abs(vdispl_max)]);
 min_plotsize = min(plot_xmax - plot_xmin,plot_ymax - plot_ymin);
 deform_fac = (min_plotsize / max_displ) * 0.15;

 % Create figure for mesh and deformed mesh plot and get its handle.
 % Locate figure at the up left corner of screen.
 fig_deform = figure;
 fig_deform_pos = get( fig_deform, 'Position' );
 fig_deform_pos(1) = 0;
 set( fig_deform, 'Position', fig_deform_pos );
 title_text = sprintf( 'Mesh and deformed mesh. Deformed factor: %s', ...
                       num2str(deform_fac) );
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
 fig_strbar = figure;
 fig_strbar_pos = get( fig_strbar, 'Position' );
 fig_strbar_pos(1) = screen_sizes(3) - fig_strbar_pos(3);
 set( fig_strbar, 'Position', fig_strbar_pos );
 title( 'Principal stress directions' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 hold on;

 % Create figure for sigma x plot and get its handle.
 % Locate figure at the second level left side of screen.
 fig_sx = figure;
 fig_sx_pos = get( fig_sx, 'Position' );
 fig_sx_pos(1) = 0;
 fig_sx_pos(2) = (screen_sizes(4) - fig_sx_pos(4))/2;
 set( fig_sx, 'Position', fig_sx_pos );
 title( 'Sigma X stress component' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 caxis([sx_nodeextrap_min sx_nodeextrap_max]);
 colorbar;
 hold on;

 % Create figure for sigma y plot and get its handle.
 % Locate figure at the second level right side of screen.
 fig_sy = figure;
 fig_sy_pos = get( fig_sy, 'Position' );
 fig_sy_pos(1) = screen_sizes(3) - fig_sy_pos(3);
 fig_sy_pos(2) = (screen_sizes(4) - fig_sy_pos(4))/2;
 set( fig_sy, 'Position', fig_sy_pos );
 title( 'Sigma Y stress component' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 caxis([sy_nodeextrap_min sy_nodeextrap_max]);
 colorbar;
 hold on;

 % Create figure for tau xy plot and get its handle.
 % Locate figure at the third level left side of screen.
 fig_txy = figure;
 fig_txy_pos = get( fig_txy, 'Position' );
 fig_txy_pos(1) = 0;
 fig_txy_pos(2) = (screen_sizes(4) - fig_txy_pos(4))/4;
 set( fig_txy, 'Position', fig_txy_pos );
 title( 'Tau XY stress component' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 caxis([txy_nodeextrap_min txy_nodeextrap_max]);
 colorbar;
 hold on;

 % Create figure for sigma 1 plot and get its handle.
 % Locate figure at the third level right side of screen.
 fig_s1 = figure;
 fig_s1_pos = get( fig_s1, 'Position' );
 fig_s1_pos(1) = screen_sizes(3) - fig_s1_pos(3);
 fig_s1_pos(2) = (screen_sizes(4) - fig_s1_pos(4))/4;
 set( fig_s1, 'Position', fig_s1_pos );
 title( 'Maximum principal stress' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 caxis([s1_nodeextrap_min s1_nodeextrap_max]);
 colorbar;
 hold on;

 % Create figure for sigma 2 plot and get its handle.
 % Locate figure at the forth level left side of screen.
 fig_s2 = figure;
 fig_s2_pos = get( fig_s2, 'Position' );
 fig_s2_pos(1) = 0;
 fig_s2_pos(2) = 0;
 set( fig_s2, 'Position', fig_s2_pos );
 title( 'Minimum principal stress' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 caxis([s2_nodeextrap_min s2_nodeextrap_max]);
 colorbar;
 hold on;

 % Create figure for tau max. plot and get its handle.
 % Locate figure at the forth level right side of screen.
 fig_tmax = figure;
 fig_tmax_pos = get( fig_tmax, 'Position' );
 fig_tmax_pos(1) = screen_sizes(3) - fig_tmax_pos(3);
 fig_tmax_pos(2) = 0;
 set( fig_tmax, 'Position', fig_tmax_pos );
 title( 'Maximum shear stress' );
 set( gca,'DataAspectRatio',[1 1 1] );
 axis([plot_xmin plot_xmax plot_ymin plot_ymax]);
 caxis([tmax_nodeextrap_min tmax_nodeextrap_max]);
 colorbar;
 hold on;

end
