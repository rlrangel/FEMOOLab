% File to include and initialize global variables

% type of analysis
global PLANE_STRESS        % plane stress analysis
global PLANE_STRAIN        % plane strain analysis
global AXISYMMETRIC        % axisymmetric analysis

PLANE_STRESS        = 0;
PLANE_STRAIN        = 1;
AXISYMMETRIC        = 2;

% type of element
global TRIA3               % linear triangular element
global QUAD4               % bilinear quadrilateral element
global TRIA6               % quadratic triangular element
global QUAD8               % serendipity quadratic quadrilateral element

TRIA3               = 0;
QUAD4               = 1;
TRIA6               = 2;
QUAD8               = 3;

% type of finite element edge (used internally)
global LINE2               % linear edge (used internally)
global LINE3               % quadratic edge (used internally)

LINE2               = 0;
LINE3               = 1;

% type of Gauss quadrature: line, triang., or quad.
global LINE_QUADRATURE     % line quadrature (used internally)
global TRIA_QUADRATURE     % triangular quadrature
global QUAD_QUADRATURE     % quadrilateral quadrature

LINE_QUADRATURE     = 0;
TRIA_QUADRATURE     = 1;
QUAD_QUADRATURE     = 2;

% type of response
global response_type
global DX_CONTOUR          % displac. x contour
global DY_CONTOUR          % displac. y contour
global SX_CONTOUR          % sigma x contour
global SY_CONTOUR          % sigma y contour
global TXY_CONTOUR         % tau xy contour
global S1_CONTOUR          % sigma 1 contour
global S2_CONTOUR          % sigma 2 contour
global TMAX_CONTOUR        % tau max. contour

 DX_CONTOUR          = 0;
 DY_CONTOUR          = 1;
 SX_CONTOUR          = 2;
 SY_CONTOUR          = 3;
 TXY_CONTOUR         = 4;
 S1_CONTOUR          = 5;
 S2_CONTOUR          = 6;
 TMAX_CONTOUR        = 7;

% Post-processing variables and data
global fig_deform          % deformed mesh figure id
global fig_strbar          % stress bar figure handle
global fig_sx              % sigma x figure handle
global fig_sy              % sigma y figure handle
global fig_txy             % tau xy figure handle
global fig_s1              % sigma 1 figure handle
global fig_s2              % sigma 2 figure handle
global fig_tmax            % tau max. figure handle

global plot_xmin           % min. x coordinate of axes (figure)
global plot_xmax           % max. x coordinate of axes (figure)
global plot_ymin           % min. y coordinate of axes (figure)
global plot_ymax           % max. y coordinate of axes (figure)

global gauss_stress_order  % order of Gauss quadrature for stress eval.
global n_gaussstress_pts   % number of gauss points for stress eval.
global TR                  % gauss stress to node stress transf. mtx.
global node_adjelems       % number of adjacent elems. of each node
 
global udispl              % vector of nodal displacements in x dir.
global udispl_min          % min. value of nodal displacements in x dir.
global udispl_max          % max. value of nodal displacements in x dir.
global vdispl              % vector of nodal displacements in y dir.
global vdispl_min          % min. value of nodal displacements in y dir.
global vdispl_max          % max. value of nodal displacements in y dir.
global deform_fac          % factor for drawing deformed mesh;

global x_gp                % vector of gauss point x coordinates
global y_gp                % vector of gauss point y coordinates
global sx_gp               % sigma x gauss points stress array
global sx_gp_min           % sigma x gauss points min. stress
global sx_gp_max           % sigma x gauss points max. stress
global sy_gp               % sigma y gauss points stress array
global sy_gp_min           % sigma y gauss points min. stress
global sy_gp_max           % sigma y gauss points max. stress
global txy_gp              % tau xy gauss points stress array
global txy_gp_min          % tau xy gauss points min. stress
global txy_gp_max          % tau xy gauss points max. stress
global s1_gp               % sigma 1 gauss points stress array
global s1_gp_min           % sigma 1 gauss points min. stress
global s1_gp_max           % sigma 1 gauss points max. stress
global s2_gp               % sigma 2 gauss points stress array
global s2_gp_min           % sigma 2 gauss points min. stress
global s2_gp_max           % sigma 2 gauss points max. stress
global tmax_gp             % tau max. gauss points stress array
global tmax_gp_min         % tau max. gauss points min. stress
global tmax_gp_max         % tau max. gauss points max. stress
global s1x_gp              % vector of s1 vector gauss x components
global s1y_gp              % vector of s1 vector gauss y components
global s2x_gp              % vector of s2 vector gauss x components
global s2y_gp              % vector of s2 vector gauss y components

global sx_elemextrap       % sigma x element node extrap. stress array
global sx_elemextrap_min   % sigma x element node extrap. min. stress
global sx_elemextrap_max   % sigma x element node extrap. max. stress
global sy_elemextrap       % sigma y element node extrap. stress array
global sy_elemextrap_min   % sigma y element node extrap. min. stress
global sy_elemextrap_max   % sigma y element node extrap. max. stress
global txy_elemextrap      % tau xy element node extrap. stress array
global txy_elemextrap_min  % tau xy element node extrap. min. stress
global txy_elemextrap_max  % tau xy element node extrap. max. stress
global s1_elemextrap       % sigma 1 element node extrap. stress array
global s1_elemextrap_min   % sigma 1 element node extrap. min. stress
global s1_elemextrap_max   % sigma 1 element node extrap. max. stress
global s2_elemextrap       % sigma 2 element node extrap. stress array
global s2_elemextrap_min   % sigma 2 element node extrap. min. stress
global s2_elemextrap_max   % sigma 2 element node extrap. max. stress
global tmax_elemextrap     % tau max. element node extrap. stress array
global tmax_elemextrap_min % tau max. element node extrap. min stress
global tmax_elemextrap_max % tau max. element node extrap. max stress

global sx_nodeextrap       % sigma x extrap. node smoothed stress array
global sx_nodeextrap_min   % sigma x extrap. node smoothed min. stress
global sx_nodeextrap_max   % sigma x extrap. node smoothed max. stress
global sy_nodeextrap       % sigma y extrap. node smoothed stress array
global sy_nodeextrap_min   % sigma y extrap. node smoothed min. stress
global sy_nodeextrap_max   % sigma y extrap. node smoothed max. stress
global txy_nodeextrap      % tau xy extrap. node smoothed stress array
global txy_nodeextrap_min  % tau xy extrap. node smoothed min. stress
global txy_nodeextrap_max  % tau xy extrap. node smoothed max. stress
global s1_nodeextrap       % sigma 1 extrap. node smoothed stress array
global s1_nodeextrap_min   % sigma 1 extrap. node smoothed min. stress
global s1_nodeextrap_max   % sigma 1 extrap. node smoothed max. stress
global s2_nodeextrap       % sigma 2 extrap. node smoothed stress array
global s2_nodeextrap_min   % sigma 2 extrap. node smoothed min. stress
global s2_nodeextrap_max   % sigma 2 extrap. node smoothed max. stress
global tmax_nodeextrap     % tau max. extrap. node smoothed stress array
global tmax_nodeextrap_min % tau max. extrap. node smoothed min. stress
global tmax_nodeextrap_max % tau max. extrap. node smoothed max. stress
