% File to include global variables

% General types
global analysis_type       % type of analysis
global PLANE_STRESS
global PLANE_STRAIN
global AXISYMMETRIC

global element_type        % type of finite element
global TRIA3               %  linear triangular element
global QUAD4               %  bilinear quadrilateral element
global TRIA6               %  quadratic triangular element
global QUAD8               %  serendipity quadratic quadrilateral element

global edge_type           % type of finite element edge (used internally)
global LINE2               %  linear edge
global LINE3               %  quadratic edge

global gauss_type          % type of Gauss quadrature: line, triang., or quad.
global LINE_QUADRATURE     %  line quadrature (used internally)
global TRIA_QUADRATURE     %  triangular quadrature
global QUAD_QUADRATURE     %  quadrilateral quadrature

% Global size parameters
global ndof                % number of degrees-of-freedom per node
global nnp                 % number of nodal points
global nel                 % number of elements
global nen                 % number of element nodes
global nedgen              % number of edge nodes
global neq                 % number of equations
global neqfree             % number of equations of free d.o.f.
global neqfixed            % number of equations of fixed d.o.f.

% Element properties variables and arrays
global gauss_order         % order of Gauss quadrature for stiffness matrix
global nmat                % number of materials
global Elasticity          % array of material elasticity coeff. values
global Poisson             % array of material poisso ratio values
global emat                % array of element material id
global nthk                % number of thicknesses
global Thickness           % array of thickness values
global ethk                % array of element thickness id

% Element nodal connectivity array and global d.o.f. numbering matrix
global IEN                 % Element nodal connectivity array: IEN(nen,nel)
                           %  stores global node number for
                           %  each local node of each element
global ID                  % global d.o.f. numbering matrix:
                           %  ID(ndof,nnp)
                           %  stores reordered global d.o.f. number for
                           %  each node
global gle                 % element gather vector (auxiliary variable):
                           %  stores eqn. numbers for element d.o.f.'s

% Nodal coordinates
global x                   % array of nodal x coordinates
global y                   % array of nodal y coordinates

% Essential Boundary Conditions (B.C.) information (support conditions)
global nnebc               % number of nodes with essential B.C.
global ebc_node            % array of nodes with essential B.C.
global ebc_flag            % array of nodal B.C. flags:
                           %  ebc_flag(ndof,nnebc)
                           %   ebc_flag = 0 --> free d.o.f.
                           %   ebc_flag = 1 --> fixed d.o.f
global nnprescdispl        % number of nodes with prescribed displac.
global prescdispl_node     % array of nodes with prescribed displac.
global prescdispl_value    % array of prescribed displac. values:
                           %  prescdispl_val(ndof,nnprescdispl)
                           %   Dx = prescdispl_value(1,j), j=1:nnprescdispl
                           %   Dy = prescdispl_value(2,j)

% Arrays of point (nodal), edge (element side) and area (element) loads
global n_pointload         % number of point loads
global pointload_node      % array of nodes with point load
                           %  pointload_node(n_pointload)
                           %   node = pointload_node(j), j=1:n_pointload
global pointload_value     % array of point load values:
                           %  pointload_value(ndof,n_pointload)
                           %   Fx = pointload_value(1,j)
                           %   Fy = pointload_value(2,j)
global n_edgeload          % number of edge uniform loads
global edgeload_side       % array of element sides with edge load
                           %  edgeload_side(3,n_edgeload)
                           %   elem = edgeload_side(1,j), j=1:n_edgeload
                           %   nod1 = edgeload_side(2,j), init side node
                           %   nod2 = edgeload_side(3,j), end side node
global edgeload_value      % array of edge uniform load values
                           %  edgeload_value(ndof,n_edgeload)
                           %   px = edgeload_value(1,j)
                           %   py = edgeload_value(2,j)
global n_areaload          % number of area uniform loads
global areaload_elem       % array of elements with area uniform load
                           %  area_load_elem(n_areaload)
                           %   elem = areaload_elem(j), j=1:n_areaload
global areaload_value      % array of elements with area load
                           %  areaload_value(ndof,n_areaload)
                           %   qx = areaload_value(1,j)
                           %   qy = areaload_value(2,j)

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

global response_type       % type of response
global DX_CONTOUR
global DY_CONTOUR
global SX_CONTOUR
global SY_CONTOUR
global TXY_CONTOUR
global S1_CONTOUR
global S2_CONTOUR
global TMAX_CONTOUR

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

