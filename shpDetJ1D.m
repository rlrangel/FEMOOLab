function detJ = shpDetJ1D(r,C)
% Evaluates 1D (edge) Jacobian determinant.
% Evaluates determinant of Jacobian matrix of parametric to cartesian
% transformation at point (r) in parametric coordinate.
% Input arguments:
%  r:    parametric coordinate value
%  C:    nodal coordinate matrix
% Output arguments:
%  detJ: determinant of Jacobian matrix evaluated at given position

 include_gblrefs;

 % Calculate GradN matrix: derivative of shape functions in parametric
 % coordinate.

 GradN = zeros(1,nedgen);

 if(edge_type == LINE2)     % linear 1D edge
 % 
 %                          (parametric coord.)
 %                      -----> r
 %                     |
 %           1 +---------------+ 2
 %
  GradN(1,1) = -0.5;
  GradN(1,2) =  0.5;

 elseif(edge_type == LINE3)     % quadratic 1D edge
 % 
 %                          (parametric coord.)
 %                     -----> r
 %                     |
 %           1 +-------+-------+ 2
 %                     3
 %
  GradN(1,1) = -0.5 + r;
  GradN(1,2) =  0.5 + r;
  GradN(1,3) = -2.0*r;
 end

 % Compute Jacobian matrix and its determinant
 J    = GradN*C;
 detJ = sqrt(J(1)*J(1) + J(2)*J(2));

end
