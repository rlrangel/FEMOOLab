function ke = elemStiffMtx(e)
% Generate stiffness matrix for a given element.
% Input arguments:
%  e:   index of target element

 include_gblrefs;

 IENe = IEN(:,e);                   % extract local connectivity information
 C = [x(IENe); y(IENe)]';           % extract element nodes (x,y) coords.
 [ngp,w,gp] = gauss(gauss_type,gauss_order); % get gauss points and weights 

 ke = zeros(nen*ndof,nen*ndof);     % init element stiffness matrix

 thk = Thickness(ethk(e));          % get element thickness 
 E = anmEMtx(e);                    % get material constituive matrix

 for i = 1:ngp                      % loop over gauss integration points

  r = gp(1,i);                      % get gauss point parametric
  s = gp(2,i);                      % coordinates

  [B,detJ] = anmBMtx(r,s,C);        % get B matrix and determinant of
                                    % jacobian evaluated at this gauss point

  ke = ke + w(i)*detJ*(B'*thk*E*B); % accumulate gauss point contribution
                                    % to element stiffness matrix

 end

end
