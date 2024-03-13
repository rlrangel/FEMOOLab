function input_patchtest_q4
% Input data for patch test with 4 QUAD4 elements in a irregular mesh.
% There are prescribed displacaments in x direction and
% applied tractions in y direction.
%
%                                8
%                         7 +----+--------------+ 9
%                           |     \             |
%                           |      \     ---    |
%                           |  ---  \   | 4 |   |
%                           | | 3 |  \   ---    |
%                           |  ---    \        -+ 6
%                           |          \ 5  --/ |
%                         4 +-----------+--/    |
%                           |   ---    /        |
%                           |  | 1 |  /    ---  |
%                           |   ---  /    | 2 | |
%                           |       /      ---  |
%                         1 +------+------------+ 3
%                                  2

 include_gblrefs;

 analysis_type = PLANE_STRESS; % plane stress analysis
 element_type = QUAD4;         % bilinear quadrilateral element
 edge_type = LINE2;            % linear edge
 gauss_type = QUAD_QUADRATURE; % type of Gauss quadrature

 % Global size parameters
 ndof = 2;                    % number of degrees-of-freedom per node
 nnp = 9;                     % number of nodal points
 nel = 4;                     % number of elements
 nen = 4;                     % number of element nodes
 nedgen = 2;                  % number of edge nodes
 neq = nnp*ndof;              % number of equations

 % Element properties variables and arrays
 gauss_order = 2;             % order of Gauss quadrature for stiffness matrix
 gauss_stress_order = 1;      % order of Gauss quadrature for stress evaluation
 nmat = 1;                    % number of materials
 Elasticity = ones(1,nmat)*10000.0;% array of material elasticity coeff values
 Poisson = ones(1,nmat)*0.25; % array of material poisso ratio values
 emat = ones(1,nel);          % array of element material id
 nthk = 1;                    % number of thicknesses
 Thickness = ones(1,nthk)*1.0;% array of thickness values
 ethk = ones(1,nel);          % array of element thickness id

 % Element nodal connectivity array
 IEN = zeros(nen,nel);
 IEN =  [1   2   4   5
	 2   3   5   6
	 5   6   8   9
	 4   5   7   8];

 % Nodal coordinates
 x = zeros(1,nnp);
 y = zeros(1,nnp);
 x = [ 0.0   3.5  10.0   0.0   6.0  10.0   0.0   2.5  10.0];
 y = [ 0.0   0.0   0.0   4.0   4.0   6.0  10.0  10.0  10.0];

 % Essential Boundary Conditions (B.C.) information (support conditions)
 nnebc = 7;                   % number of nodes with essential B.C.
 ebc_node = zeros(1,nnebc);   % array of nodes with essential B.C.
 ebc_node = [1  2  3  4  6  7  9];
 ebc_flag = zeros(ndof,nnebc); % array of nodal B.C. flags
 ebc_flag = [1  0  1  1  1  1  1
	     1  1  1  0  0  0  0];
 nnprescdispl = 3;             % number of nodes with prescrib. diplac.
 prescdispl_node = zeros(1,nnprescdispl);     % array of nodes with presc.disp.
 prescdispl_node = [3  6  9];
 prescdispl_value = zeros(ndof,nnprescdispl); % array of presc.disp. values
 prescdispl_value = [0.005  0.005  0.005
                     0.000  0.000  0.000];

 % Arrays of point (nodal), edge (element side) and area (element) loads
 n_pointload = 0;             % number of point loads
 n_edgeload = 2;              % number of edge uniform loads
 edgeload_side = zeros(3,n_edgeload);% array of element sides with edge load
 edgeload_side = [3  4     % elem
		  7  8     % nod1
		  8  9];   % nod2
 edgeload_value = zeros(ndof,n_edgeload);% array of edge uniform load values
 edgeload_value = [ 0.0   0.0    % px
		   20.0  20.0];  % py
 n_areaload = 0;              % number of area uniform loads

end
