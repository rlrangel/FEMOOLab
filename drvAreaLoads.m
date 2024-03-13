function F = drvAreaLoads(F)
% Add area (element) equivalent nodal loads to global forcing vector.
% Input/output arguments:
%  F:    global forcing vector

 include_gblrefs;

 for j = 1:n_areaload
  e = areaload_elem(j);
  [nne,IENe,feq] = elemAreaENL(e,areaload_value(:,j));
  F = drvAssembleENL(F,nne,IENe,feq);
 end

end
