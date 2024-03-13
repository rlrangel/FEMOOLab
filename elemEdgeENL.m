function [nne,IENe,feq] = elemEdgeENL(e,nod1,nod2,p)
% Generate equivalent nodal load vector (feq) for an edge uniform distributed load (px,py) on a given element side.
% Input arguments:
%  e:      element number
%  nod1:   first node on element edge
%  nod2:   last node on element edge
%  p:      edge uniform distributed load values: p(ndof,1)
% Output arguments:
%  nne:    number of nodes of equivalent nodal load vector
%  IENe:   local nodal connectivity array of equivalent nodal load
%  feq:    equivalent nodal load vector

 include_gblrefs;

 [nne,IENe] = shpGetEdge(e,nod1,nod2); % get nodes on element side

 C = [x(IENe); y(IENe)]';           % extract node (x,y) coords.
 gtype = LINE_QUADRATURE;           % unidimension gauss integration
 [ngp,w,gp] = gauss(gtype,gauss_order);  % get gauss points and weights 

 feq = zeros(nne*ndof,1);           % init element equiv. nodal load vector

 for i = 1:ngp                      % loop over gauss integration points

  r = gp(1,i);                      % get gauss point parametric coordinate

  N = shpNmatrix1D(r);              % get shape functions matrix
                                    % evaluated at this gauss point
  detJ = shpDetJ1D(r,C);            % get determinant of jacobian 
                                    % evaluated at this gauss point
  il = 0;
  for j = 1:nne                     % accumulate gauss point contribution
   for k = 1:ndof
    il = il + 1;
    feq(il) = feq(il) + w(i)*detJ*N(j)*p(k);
   end
  end

 end

end
