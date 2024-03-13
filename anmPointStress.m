function str = anmPointStress(r,s,C,E,d)
% Get stress components (sigma x, sigma y and tau xy) at a given point on an element.
% For axisymmetric analysis the sigma z stress component is desregarded.
% Input arguments:
%  r,s:  parametric coordinate values of point location
%  C:    nodal coordinate matrix of target element
%  E:    constituive matrix of target element
%  d:    generalized displacements/rotations for all d.o.f.'s of element
% Output arguments:
%  str:  stress components (sx,sy,txy) at target point

 include_gblrefs;

 str = zeros(3,1);                  % init stress component vector

 [B,detJ] = anmBMtx(r,s,C);         % get B matrix and determinant of
                                    % jacobian evaluated at this gauss point

 str_raw = E*B*d;                   % compute point stress components

 if(analysis_type == AXISYMMETRIC)  % skip tangential stress component

  str(1) = str_raw(1);
  str(2) = str_raw(3);
  str(3) = str_raw(4);

 else % in plane stress or plane strain raw stress vector is the target one

  str = str_raw;

 end

end
