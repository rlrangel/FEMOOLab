function posElemStressesExtrap
% Create arrays and compute node extrapolated stress components and principal stressesfor all elements.
% The nodal stress components are computed by extrapolation of gauss point 
% stress components using the TR matrix (which is a global variable that is
% assumed already created).
%
% All stress arrays are declared as global variables
% (see file include_gblrefs.m).

 include_gblrefs;

 sx_elemextrap   = zeros(nen,nel);
 sy_elemextrap   = zeros(nen,nel);
 txy_elemextrap  = zeros(nen,nel);
 s1_elemextrap   = zeros(nen,nel);
 s2_elemextrap   = zeros(nen,nel);
 tmax_elemextrap = zeros(nen,nel);

 for e = 1:nel

  sx_elemextrap(:,e)   = TR * sx_gp(:,e);
  sy_elemextrap(:,e)   = TR * sy_gp(:,e);
  txy_elemextrap(:,e)  = TR * txy_gp(:,e);
  s1_elemextrap(:,e)   = TR * s1_gp(:,e);
  s2_elemextrap(:,e)   = TR * s2_gp(:,e);
  tmax_elemextrap(:,e) = TR * tmax_gp(:,e);

 end

 sx_elemextrap_min   = min(min(sx_elemextrap));
 sx_elemextrap_max   = max(max(sx_elemextrap));
 sy_elemextrap_min   = min(min(sy_elemextrap));
 sy_elemextrap_max   = max(max(sy_elemextrap));
 txy_elemextrap_min  = min(min(txy_elemextrap));
 txy_elemextrap_max  = max(max(txy_elemextrap));
 s1_elemextrap_min   = min(min(s1_elemextrap));
 s1_elemextrap_max   = max(max(s1_elemextrap));
 s2_elemextrap_min   = min(min(s2_elemextrap));
 s2_elemextrap_max   = max(max(s2_elemextrap));
 tmax_elemextrap_min = min(min(tmax_elemextrap));
 tmax_elemextrap_max = max(max(tmax_elemextrap));

 % Clear small value responses
 if((abs(sx_elemextrap_min) < 0.00001)&&(abs(sx_elemextrap_max) < 0.00001))
  sx_elemextrap_min = 0.0;
  sx_elemextrap_max = 0.0;
  sx_elemextrap = zeros(nen,nel);
 end
 if((abs(sy_elemextrap_min) < 0.00001)&&(abs(sy_elemextrap_max) < 0.00001))
  sy_elemextrap_min = 0.0;
  sy_elemextrap_max = 0.0;
  sy_elemextrap = zeros(nen,nel);
 end
 if((abs(txy_elemextrap_min) < 0.00001)&&(abs(txy_elemextrap_max) < 0.00001))
  txy_elemextrap_min = 0.0;
  txy_elemextrap_max = 0.0;
  txy_elemextrap = zeros(nen,nel);
 end
 if((abs(s1_elemextrap_min) < 0.00001)&&(abs(s1_elemextrap_max) < 0.00001))
  s1_elemextrap_min = 0.0;
  s1_elemextrap_max = 0.0;
  s1_elemextrap = zeros(nen,nel);
 end
 if((abs(s2_elemextrap_min) < 0.00001)&&(abs(s2_elemextrap_max) < 0.00001))
  s2_elemextrap_min = 0.0;
  s2_elemextrap_max = 0.0;
  s2_elemextrap = zeros(nen,nel);
 end
 if((abs(tmax_elemextrap_min) < 0.00001)&&(abs(tmax_elemextrap_max) < 0.00001))
  tmax_elemextrap_min = 0.0;
  tmax_elemextrap_max = 0.0;
  tmax_elemextrap = zeros(nen,nel);
 end

end
