function F = drvAssembleENL(F,nne,IENe,feq)
% Add equivalent nodal load vector to the global forcing vector.
% Assemble the given equivalent nodal load vector to any term of
% the global forcing vector, incluing the terms that correspond
% to constrained d.o.f.
% Input/output arguments:
%  F:    global forcing vector
% Input arguments:
%  nne:  number of nodes of given equivalent nodal load vector
%  IENe: array of nodes of given equivalent nodal load
%  feq:  equivalent nodal load vector

 include_gblrefs;

 % Add contribution of equivalent nodal load vector to global forcing vector
 % il is the local (element) row equation number
 % ig is the global row equation number
 il = 0;
 for j = 1:nne
  for k = 1:ndof
   il = il + 1;
   ig = ID(k,IENe(j));
   F(ig) = F(ig) + feq(il);   % assemble given equiv. nodal load vector
  end
 end

end
