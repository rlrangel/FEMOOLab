function [nne,IENe,feq] = elemAreaENL(e,q)
% Generate equivalent nodal load vector (feq) for an area uniform distributed load (qx,qy) on a given element.
% Input arguments:
%  e:      element number
%  q:      area uniform distributed load values: q(ndof,1)
% Output arguments:
%  nne:    number of nodes of equivalent nodal load vector
%  IENe:   local nodal connectivity array of equivalent nodal load
%  feq:    equivalent nodal load vector

 include_gblrefs;

 nne = nen;                         % number of nodes of equiv.load vector
				    % is equal to number of element nodes

 IENe = IEN(:,e);                   % extract local connectivity information
 C = [x(IENe); y(IENe)]';           % extract node (x,y) coords.
 [ngp,w,gp] = gauss(gauss_type,gauss_order); % get gauss points and weights 

 feq = zeros(nne*ndof,1);           % init element equiv. nodal load vector

 for i = 1:ngp                      % loop over gauss integration points

  r = gp(1,i);                      % get gauss point parametric
  s = gp(2,i);                      % coordinates

  N = shpNmatrix2D(r,s);            % get shape functions matrix
                                    % evaluated at this gauss point
  detJ = shpDetJ2D(r,s,C);          % get determinant of jacobian 
                                    % evaluated at this gauss point
  il = 0;
  for j = 1:nne                     % accumulate gauss point contribution
   for k = 1:ndof
    il = il + 1;
    feq(il) = feq(il) + w(i)*detJ*N(j)*q(k);
   end
  end

 end

end
