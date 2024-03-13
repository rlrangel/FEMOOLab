function posProcessor(D)
% Visualize mesh and results of current analysis.
% Input arguments:
%  D:   global displacement vector

 include_gblrefs;

 fprintf(1,'Preparing post-processing data...\n');

 % Compute transformation matrix to transform gauss point
 % results into element node results.
 elemTRMtx;

 % Create nodal displacement response data
 udispl = D(ID(1,:));
 vdispl = D(ID(2,:));
 udispl_min = min(udispl);
 udispl_max = max(udispl);
 vdispl_min = min(vdispl);
 vdispl_max = max(vdispl);

 % Compute gauss point stress components and principal stresses.
 % Store data in global variable tables.
 posGaussStresses(D);

 % Extrapolate gauss point results to element node results.
 % Store data in global variable tables.
 posElemStressesExtrap;

 % Average smooths element node result to global node results.
 % Store data in global variable tables.
 posNodeStressesExtrap;

 fprintf(1,'Displaying results...\n');

 % Setup bounding box for displaying results.
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
 posCreateFigs;

 % Display deformed mesh.
 figure(fig_deform);
 posPlotMesh('k');
 posPlotDeformMesh('b');

 % Display principal stress vectors.
 figure(fig_strbar);
 posPlotMesh('k');
 quiver(x_gp,y_gp,s1x_gp,s1y_gp,'r');
 quiver(x_gp,y_gp,s2x_gp,s2y_gp,'b');

 % Display stress results.
 figure(fig_sx);
 response_type = SX_CONTOUR;
 posPlotNodeContour;
 posPlotMesh('k');

 figure(fig_sy);
 response_type = SY_CONTOUR;
 posPlotNodeContour;
 posPlotMesh('k');

 figure(fig_txy);
 response_type = TXY_CONTOUR;
 posPlotNodeContour;
 posPlotMesh('k');

 figure(fig_s1);
 response_type = S1_CONTOUR;
 posPlotNodeContour;
 posPlotMesh('k');

 figure(fig_s2);
 response_type = S2_CONTOUR;
 posPlotNodeContour;
 posPlotMesh('k');

 figure(fig_tmax);
 response_type = TMAX_CONTOUR;
 posPlotNodeContour;
 posPlotMesh('k');

end
