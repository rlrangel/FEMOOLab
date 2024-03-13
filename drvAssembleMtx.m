function K = drvAssembleMtx(K,e,ke)
% Assemble global stiffness matrix.
% Input/output arguments:
%  K:    global stiffness matrix
% Input arguments:
%  e:    element number
%  ke:   element local stiffness matrix

 include_gblrefs;

 % Assemble element gather vector (stores element d.o.f. eqn. numbers)
 i = 0;
 for j = 1:nen
  for k = 1:ndof
   i = i + 1;
   gle(i) = ID(k,IEN(j,e));
  end
 end

 % Add contribution of element given by its id to global stiffness matrix
 % il is the local (element) row equation number
 % jl is the local (element) column equation number
 % ig is the global row equation number
 % jg is the global column equation number
 for il = 1:nen*ndof
  ig = gle(il);
  for jl = 1:nen*ndof
   jg = gle(jl);
   K(ig,jg) = K(ig,jg) + ke(il,jl);      % assemble global stiffness matrix
  end
 end

end
