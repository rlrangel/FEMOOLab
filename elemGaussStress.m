function [ngp,str,gpc] = elemGaussStress(D,e)
% Get stress components at gauss points and gauss point cartesian coordinates for a given element.
% Input arguments:
%  D:    solution vector: generalized displacements/rotations for all d.o.f.'s
%  e:    target finite element index
% Output arguments:
%  ngp:  number of gauss points for stress evaluation
%  str:  stress components (sx,sy,txy) at each gauss point
%  gpc:  gauss point cartesian coordinates array

 include_gblrefs;

 IENe   = IEN(:,e);                  % extract local connectivity information

 % assemble element gather vector (stores element d.o.f. eqn. numbers)
 i = 0;
 for j = 1:nen
  for k = 1:ndof
   i = i + 1;
   gle(i) = ID(k,IENe(j));
  end
 end

 d = D(gle);                         % extract solution at element nodes

 C = [x(IENe); y(IENe)]';            % extract element nodes (x,y) coords.
 [ngp,w,gp] = gauss(gauss_type,gauss_stress_order); % get Gauss points and weights

 str = zeros(3,ngp);                 % init element stress component matrix
 gpc = zeros(2,ngp);                 % init element gauss point coordinates

 E = anmEMtx(e);                     % get material constituive matrix

 % compute stress components at gauss points and 
 % gauss point cartesian coordinates
 for i = 1:ngp

  r = gp(1,i);                       % get gauss point parametric
  s = gp(2,i);                       % coordinates

  N = shpNmatrix2D(r,s);             % get shape functions matrix
                                     % evaluated at this gauss point

  str(:,i) = anmPointStress(r,s,C,E,d); % compute gauss point stress components
  gpc(:,i) = N*C;                    % and gauss points cartesian coordinates  

 end

end
