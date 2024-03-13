function posNodeStressesExtrap
% Create arrays and compute extrapolated node smoothed stress components and principal stresses for all nodes.
% The nodal stress components are computed by averaging values of element 
% extrapolated nodal stress components of all elements adjacent to each node.
%
% All stress arrays are declared as global variables
% (see file include_gblrefs.m).

 include_gblrefs;

 % assemble vector with number of adjacent elements of each node
 if isempty(node_adjelems)
  node_adjelems = zeros(nnp,1);
  for e = 1:nel
   for i = 1:nen
    n = IEN(i,e);
    node_adjelems(n) = node_adjelems(n) + 1;
   end
  end
 end

 sx_nodeextrap   = zeros(nnp,1);
 sy_nodeextrap   = zeros(nnp,1);
 txy_nodeextrap  = zeros(nnp,1);
 s1_nodeextrap   = zeros(nnp,1);
 s2_nodeextrap   = zeros(nnp,1);
 tmax_nodeextrap = zeros(nnp,1);

 for e = 1:nel
  for i = 1:nen
   n = IEN(i,e);
   sx_nodeextrap(n)   = sx_nodeextrap(n)   + sx_elemextrap(i,e);
   sy_nodeextrap(n)   = sy_nodeextrap(n)   + sy_elemextrap(i,e);
   txy_nodeextrap(n)  = txy_nodeextrap(n)  + txy_elemextrap(i,e);
   s1_nodeextrap(n)   = s1_nodeextrap(n)   + s1_elemextrap(i,e);
   s2_nodeextrap(n)   = s2_nodeextrap(n)   + s2_elemextrap(i,e);
   tmax_nodeextrap(n) = tmax_nodeextrap(n) + tmax_elemextrap(i,e);
  end
 end

 for n = 1:nnp
  sx_nodeextrap(n)   = sx_nodeextrap(n)   / node_adjelems(n);
  sy_nodeextrap(n)   = sy_nodeextrap(n)   / node_adjelems(n);
  txy_nodeextrap(n)  = txy_nodeextrap(n)  / node_adjelems(n);
  s1_nodeextrap(n)   = s1_nodeextrap(n)   / node_adjelems(n);
  s2_nodeextrap(n)   = s2_nodeextrap(n)   / node_adjelems(n);
  tmax_nodeextrap(n) = tmax_nodeextrap(n) / node_adjelems(n);
 end

 sx_nodeextrap_min   = min(min(sx_nodeextrap));
 sx_nodeextrap_max   = max(max(sx_nodeextrap));
 sy_nodeextrap_min   = min(min(sy_nodeextrap));
 sy_nodeextrap_max   = max(max(sy_nodeextrap));
 txy_nodeextrap_min  = min(min(txy_nodeextrap));
 txy_nodeextrap_max  = max(max(txy_nodeextrap));
 s1_nodeextrap_min   = min(min(s1_nodeextrap));
 s1_nodeextrap_max   = max(max(s1_nodeextrap));
 s2_nodeextrap_min   = min(min(s2_nodeextrap));
 s2_nodeextrap_max   = max(max(s2_nodeextrap));
 tmax_nodeextrap_min = min(min(tmax_nodeextrap));
 tmax_nodeextrap_max = max(max(tmax_nodeextrap));

end
