function F = drvPointLoads(F)
% Add point (nodal) loads to global forcing vector.
% Add nodal load components to any term of the global forcing vector,
% incluing the terms that correspond to constrained d.o.f.
% Input/output arguments:
%  F:    global forcing vector

 include_gblrefs;

 for j = 1:n_pointload
  n = pointload_node(j);
  for k = 1:ndof
   i = ID(k,n);
   F(i) = F(i) + pointload_value(k,j);  % add nodal load component
  end
 end

end
