function [B,detJ] = anmBMtx(r,s,C)
% Assemble 2D B matrix for an element.
% Evaluates derivatives of shape functions (in physical Cartesian coordinates) 
% at point (r,s) in parametric coordinates and assembles them in the B matrix.
% In case of axisymmetric analysis, the shape function itself is also used.
% Input arguments:
%  r,s:  parametric coordinate values
%  C:    nodal coordinate matrix
% Output arguments:
%  B:    B matrix evaluated at given position

 include_gblrefs;

 % Calculate GradNpar matrix: derivative of shape functions in parametric
 % coordinates.

 GradNpar = shpGradNmatrix2D(r,s);

 % Compute Jacobian matrix and its determinant
 J    = GradNpar*C;
 detJ = det(J);

 % Calculate GradNcar matrix: derivative of shape functions in cartesian
 % coordinates.
 GradNcar = J\GradNpar;

 % Assemble B matrix
 if(analysis_type == AXISYMMETRIC)

  N = shpNmatrix2D(r,s);

  p = N*C;        % location of evaluation point
  radius = p(1);  % x coordinate is the radius in axisymmetry

  B = zeros(4,2*nen);

  for i = 1:nen
   B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0.0;
   B(2,2*i-1) = N(i)/radius;     B(2,2*i) = 0.0;
   B(3,2*i-1) = 0.0;             B(3,2*i) = GradNcar(2,i);
   B(4,2*i-1) = GradNcar(2,i);   B(4,2*i) = GradNcar(1,i);
  end

 else % plane stress and plain strain do not need shape function

  B = zeros(3,2*nen);
  for i = 1:nen
   B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0.0;
   B(2,2*i-1) = 0.0;             B(2,2*i) = GradNcar(2,i);
   B(3,2*i-1) = GradNcar(2,i);   B(3,2*i) = GradNcar(1,i);
  end

 end

end
