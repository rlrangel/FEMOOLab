function [K,F,D] = preProcessor
% Preprocessing input data and sets up mesh information.
% Output arguments:
%  K:   global stiffness matrix (returns initialized with null values)
%  F:   global forcing vector (returns initialized with null values)
%  D:   global displacement vector (returns initialized with null values
%       and with known support settlement values)

 include_gblrefs;

 % input file to include all variables

 % get neutral file full file name
 filterspec = {'*.nf';'*.dat';'*.pos'};
 [filename,pathname] = uigetfile(filterspec,'Elasticity2D - Input file');
 if isequal(filename,0)
  input_patchtest_q4;
 else
  fullname = strcat(pathname,filename);
  preReadNF(fullname);
 end

 fprintf(1,'Pre-processing...\n');

 % Initialize global matrix and vectors
 K = zeros(neq,neq);          % initialize global stiffness matrix
 F = zeros(neq,1);            % initialize global system forcing vector
 D = zeros(neq,1);            % initialize global system displacement vector

 % Initialize auxiliary variables
 gle = zeros(nen*ndof,1);     % element gather vector:
			      %  stores eqn. numbers for element d.o.f.'s

 % generate global d.o.f numbering matrix
 nodeSetupDOFNum;

 % store prescribed displacements in global displacement vector
 % (known support settlement values)
 for i = 1:nnprescdispl
  n = prescdispl_node(i);  % get the node id
  for k = 1:ndof
   if(ID(k,n) > neqfree)   % check to see whether this d.o.f is really fixed
    D(ID(k,n)) = prescdispl_value(k,i);
  end
 end

end
