function init_constants
% Initialize constant global variables.

 include_gblrefs;

 PLANE_STRESS        = 0;    % plane stress analysis
 PLANE_STRAIN        = 1;    % plane strain analysis
 AXISYMMETRIC        = 2;    % axisymmetric analysis

 TRIA3               = 0;    % linear triangular element
 QUAD4               = 1;    % bilinear quadrilateral element
 TRIA6               = 2;    % quadratic triangular element
 QUAD8               = 3;    % serendipity quadratic quadrilateral element

 LINE2               = 0;    % linear edge (used internally)
 LINE3               = 1;    % quadratic edge (used internally)

 LINE_QUADRATURE     = 0;    % line quadrature (used internally)
 TRIA_QUADRATURE     = 1;    % triangular quadrature
 QUAD_QUADRATURE     = 2;    % quadrilateral quadrature

 DX_CONTOUR          = 0;    % displac. x contour
 DY_CONTOUR          = 1;    % displac. y contour
 SX_CONTOUR          = 2;    % sigma x contour
 SY_CONTOUR          = 3;    % sigma y contour
 TXY_CONTOUR         = 4;    % tau xy contour
 S1_CONTOUR          = 5;    % sigma 1 contour
 S2_CONTOUR          = 6;    % sigma 2 contour
 TMAX_CONTOUR        = 7;    % tau max. contour

end
