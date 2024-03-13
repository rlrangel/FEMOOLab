function posGaussStresses(D)
% Create arrays and compute gauss stress components and principal stresses for all elements.
% Input arguments:
%  D:    solution vector: generalized displacements/rotations for all d.o.f.'s
%
% All gauss point location and stress arrays are declared as global variables
% (see file incinclude_gblrefs.m).

 include_gblrefs;

 x_gp      = zeros(n_gaussstress_pts*nel,1);
 y_gp      = zeros(n_gaussstress_pts*nel,1);
 sx_gp     = zeros(n_gaussstress_pts,nel);
 sy_gp     = zeros(n_gaussstress_pts,nel);
 txy_gp    = zeros(n_gaussstress_pts,nel);
 s1_gp     = zeros(n_gaussstress_pts,nel);
 s2_gp     = zeros(n_gaussstress_pts,nel);
 tmax_gp   = zeros(n_gaussstress_pts,nel);
 s1x_gp    = zeros(n_gaussstress_pts*nel,1);
 s1y_gp    = zeros(n_gaussstress_pts*nel,1);
 s2x_gp    = zeros(n_gaussstress_pts*nel,1);
 s2y_gp    = zeros(n_gaussstress_pts*nel,1);

 npts = 0;
 for e = 1:nel

  [ngp,str,gpc] = elemGaussStress(D,e);

  sx_gp(:,e)   = str(1,:);
  sy_gp(:,e)   = str(2,:);
  txy_gp(:,e)  = str(3,:);

  for i = 1:ngp
   npts = npts + 1;
   x_gp(npts) = gpc(1,i);
   y_gp(npts) = gpc(2,i);
   [prc,thetap] = posPrincStress(str(:,i));
   s1_gp(i,e)   = prc(1);
   s2_gp(i,e)   = prc(2);
   tmax_gp(i,e) = prc(3);
   s1x_gp(npts) = s1_gp(i,e)*cos(thetap);
   s1y_gp(npts) = s1_gp(i,e)*sin(thetap);
   s2x_gp(npts) = s2_gp(i,e)*cos(thetap+(pi/2.0));
   s2y_gp(npts) = s2_gp(i,e)*sin(thetap+(pi/2.0));
  end

 end

 sx_gp_min   = min(min(sx_gp));
 sx_gp_max   = max(max(sx_gp));
 sy_gp_min   = min(min(sy_gp));
 sy_gp_max   = max(max(sy_gp));
 txy_gp_min  = min(min(txy_gp));
 txy_gp_max  = max(max(txy_gp));
 s1_gp_min   = min(min(s1_gp));
 s1_gp_max   = max(max(s1_gp));
 s2_gp_min   = min(min(s2_gp));
 s2_gp_max   = max(max(s2_gp));
 tmax_gp_min = min(min(tmax_gp));
 tmax_gp_max = max(max(tmax_gp));

 % Clear small value responses
 if((abs(sx_gp_min) < 0.00001)&&(abs(sx_gp_max) < 0.00001))
  sx_gp_min = 0.0;
  sx_gp_max = 0.0;
  sx_gp = zeros(n_gaussstress_pts,nel);
 end
 if((abs(sy_gp_min) < 0.00001)&&(abs(sy_gp_max) < 0.00001))
  sy_gp_min = 0.0;
  sy_gp_max = 0.0;
  sy_gp = zeros(n_gaussstress_pts,nel);
 end
 if((abs(txy_gp_min) < 0.00001)&&(abs(txy_gp_max) < 0.00001))
  txy_gp_min = 0.0;
  txy_gp_max = 0.0;
  txy_gp = zeros(n_gaussstress_pts,nel);
 end
 if((abs(s1_gp_min) < 0.00001)&&(abs(s1_gp_max) < 0.00001))
  s1_gp_min = 0.0;
  s1_gp_max = 0.0;
  s1_gp = zeros(n_gaussstress_pts,nel);
 end
 if((abs(s2_gp_min) < 0.00001)&&(abs(s2_gp_max) < 0.00001))
  s2_gp_min = 0.0;
  s2_gp_max = 0.0;
  s2_gp = zeros(n_gaussstress_pts,nel);
 end
 if((abs(tmax_gp_min) < 0.00001)&&(abs(tmax_gp_max) < 0.00001))
  tmax_gp_min = 0.0;
  tmax_gp_max = 0.0;
  tmax_gp = zeros(n_gaussstress_pts,nel);
 end

end
