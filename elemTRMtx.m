function elemTRMtx
% Compute gauss-to-node results transfer matrix (TR matrix).
% Refs.:
% Hinton & Campbell, "Local and Global Smoothing of Discontinous Finite
% Element Functions using a Least Squares Method", Int. J. Num. Meth. Engng.,
% Vol. 8, pp. 461-480, 1974.
% Burnett, D.S., "Finite Element Analysis - From Concepts to Applications",
% Addison-Wesley, 1987.
% Martha, L.F., "Notas de Aula do Curso CIV 2118 - Metodo dos Elementos
% Finitos", 1994.
% This functions loads global variables TR and n_gaussstress_pts
% (see file include_gblrefs.m).

 include_gblrefs;

 % Get gauss points parametric coordinates for stress evaluation quadratura
 [ngp,w,gp] = gauss(gauss_type,gauss_stress_order);

 if(gauss_stress_order == 1)

  TR = ones(nen,1);

 else

  % Compute S matrix that defines the coefficients of the smoothing stress
  % plane that fits the gauss point stress values in a least square sence.

  P = zeros(3,3);
  P(1,1) = ngp;
  P(1,2) = sum(gp(1,:));
  P(1,3) = sum(gp(2,:));
  P(2,1) = P(1,2);
  P(2,2) = gp(1,:) * gp(1,:)';
  P(2,3) = gp(1,:) * gp(2,:)';
  P(3,1) = P(1,3);
  P(3,2) = P(2,3);
  P(3,3) = gp(2,:) * gp(2,:)';

  Q = ones(3,ngp);
  Q(2,:) = gp(1,:);
  Q(3,:) = gp(2,:);

  S = inv(P) * Q;

  % Compute the nodal stress evaluation matrix, which is obtained 
  % using the nodal parametric coordinates in the smoothing stress
  % plane equation.

  pt = shpNodeParCoords2D;   % Get element node parametric coordinates

  E = ones(nen,3);
  E(:,2) = pt(1,:)';
  E(:,3) = pt(2,:)';

  % Compute the TR matrix and save number of gauss points for stress eval.

  TR = E * S;

 end

 n_gaussstress_pts = ngp;

end
