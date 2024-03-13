function F = drvEdgeLoads(F)
% Add edge (element side) equivalent nodal loads to global forcing vector.
% Input/output arguments:
%  F:    global forcing vector

 include_gblrefs;

 for j = 1:n_edgeload
  e = edgeload_side(1,j);
  nod1 = edgeload_side(2,j);
  nod2 = edgeload_side(3,j);
  [nne,IENe,feq] = elemEdgeENL(e,nod1,nod2,edgeload_value(:,j));
  F = drvAssembleENL(F,nne,IENe,feq);
 end

end
