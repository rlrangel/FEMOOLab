function nodeSetupDOFNum
% Setup global d.o.f numbering matrix.
% Free d.o.f.'s are numbered first.
% Counts total number of equations of free d.o.f.'s and total number
% of equations of fixed d.o.f.'s, and stores in global variables.

 include_gblrefs;

 ID = zeros(ndof,nnp);

 % Compute total number of fixed d.o.f.
 % Initialize ID matrix:
 %  if ID(k,n) = 0, d.o.f. k of node n is free
 %  if ID(k,n) = 1, d.o.f. k of node n is fixed
 neqfixed = 0;
 for i = 1:nnebc
  n = ebc_node(i);
  for k = 1:ndof
   if( ebc_flag(k,i) == 1 )
    neqfixed = neqfixed + 1;
    ID(k,n) = 1;
   end
  end
 end

 % Compute total number of free d.o.f.
 neqfree = neq - neqfixed;

 % Assemble ID matrix:
 %  ID(k,n) = equation number of d.o.f. k of node n
 %  Free d.o.f.'s have the initial numbering
 countl = 0;
 countf = neqfree;
 for n = 1:nnp
  for k = 1:ndof
   if(ID(k,n) == 0)
    countl = countl + 1;
    ID(k,n) = countl;
   else
    countf = countf + 1;
    ID(k,n) = countf;
   end
  end
 end

end
