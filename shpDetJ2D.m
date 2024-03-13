function detJ = shpDetJ2D(r,s,C)
% Evaluates 2D Jacobian determinant.
% Evaluates determinant of Jacobian matrix of parametric to cartesian
% transformation at point (r,s) in parametric coordinates.
% Input arguments:
%  r,s:  parametric coordinate values
%  C:    nodal coordinate matrix
% Output arguments:
%  detJ: determinant of Jacobian matrix evaluated at given position

 include_gblrefs;

 % Calculate GradNpar matrix: derivative of shape functions in parametric
 % coordinates.

 GradNpar = shpGradNmatrix2D(r,s);

 % Compute Jacobian matrix and its determinant
 J    = GradNpar*C;
 detJ = det(J);

end
