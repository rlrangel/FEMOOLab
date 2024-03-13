function N = shpNmatrix1D(r)
% Evaluates 1D shape function matrix.
% Evaluates shape functions (in parametric coordinate) at point (r)
% for 1D elements.
% Input arguments:
%  r:    parametric coordinate value
% Output arguments:
%  N:    shape function matrix evaluated at given position

 include_gblrefs;

 N = zeros(1,nedgen);

 if(edge_type == LINE2)         % linear 1D edge
 % 
 %                          (parametric coord.)
 %                      -----> r
 %                     |
 %           1 +---------------+ 2
 % 
  N(1) = 0.5*(1.0-r);
  N(2) = 0.5*(1.0+r);

 else if(edge_type == LINE3)    % quadratic 1D edge
 % 
 %                          (parametric coord.)
 %                     -----> r
 %                     |
 %           1 +-------+-------+ 2
 %                     3
 % 
  N(3) = 1.0-r*r;
  N(1) = 0.5*(1.0-r) - 0.5*N(3);
  N(2) = 0.5*(1.0+r) - 0.5*N(3);
 end

end
