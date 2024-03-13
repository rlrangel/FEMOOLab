function preReadNF(arquivo)
% This fuction reads a neutral-format file with input data for a finite-element model for this program.
% The neutral format is described in www.tecgraf.puc-rio.br/neutralfile.
% It is assumed that all nodes and elements are numbered consecutively from
% one to the total number of nodes and elements.
% It is assumed that there is only one type of element and only one
% type of gauss integration order in each finite element model.
% It is assumed that there is a single load case.
% See global variables in file include_gbl_refs.m.

 include_gblrefs;

 % Initialize some global variable with default values.
 analysis_type = PLANE_STRESS; % type of analysis: plane stress
 ndof = 2;                     % number of degrees-of-freedom per node
 nnprescdispl = 0;             % initialize number of nodes with presc.displ.
 n_pointload = 0;              % initialize number of point loads
 n_edgeload = 0;               % initialize number of edge uniform loads
 n_areaload = 0;               % initialize number of area uniform loads

 % Open neutral-format file.
 fid = fopen(arquivo,'rt');

 flag = 1;

 while flag
  % Get file line.
  tline = fgetl(fid);
  % Get rid of blank spaces.
  string = deblank(tline);
    
  % Check for a neutral-format tag string.
  switch string
   case '%HEADER.ANALYSIS'
    ReadHeaderAnalysis(fid);
   case '%NODE.COORD'
    ReadNodeCoord(fid);
   case '%NODE.SUPPORT'
    ReadNodeSupport(fid);
   case '%MATERIAL'
    ReadMaterial(fid);
   case '%MATERIAL.ISOTROPIC'
    ReadMaterialIsotropic(fid);
   case '%THICKNESS'
    ReadThickness(fid);
   case '%INTEGRATION.ORDER'
    ReadIntegrationOrder(fid);
   case '%ELEMENT'
    ReadElement(fid);
   case '%ELEMENT.T3'
    ReadElementT3(fid);
   case '%ELEMENT.Q4'
    ReadElementQ4(fid);
   case '%ELEMENT.T6'
    ReadElementT6(fid);
   case '%ELEMENT.Q8'
    ReadElementQ8(fid);
   case '%LOAD.CASE.NODAL.FORCES'
    ReadLoadCaseNodalForces(fid);
   case '%LOAD.CASE.NODAL.DISPLACEMENT'
    ReadLoadCaseNodalDisplacement(fid);
   case '%LOAD.CASE.LINE.FORCE.UNIFORM'
    ReadLoadCaseLineForceUniform(fid);
   case '%LOAD.CASE.DOMAIN.FORCE.UNIFORM'
    ReadLoadCaseDomainForceUniform(fid);
  end
    
  if strcmp(string,'%END')
   flag = 0;
  end
 end

 % Check to see whether a model thickness was specified.
 % If not, specify a single unit thickness for all elements.
 if(nthk == 0)
  nthk = 1;
  Thickness = ones(nthk);
 end

 fclose(fid);
end

%----------------------------------
function ReadHeaderAnalysis(fid)
 include_gblrefs;

 % Get file line.
 tline = fgetl(fid);
 % Get rid of blank spaces.
 string = deblank(tline);

 % Get target analysis type.
 switch string
  case '''plane_stress'''
   analysis_type = PLANE_STRESS;
  case '''plane_strain'''
   analysis_type = PLANE_STRAIN;
  case '''axisymmetric'''
   analysis_type = AXISYMMETRIC;
 end
end

%-----------------------------------
function ReadNodeCoord(fid)
 include_gblrefs;

 % Store total number of nodes of model.
 nnp = fscanf(fid,'%f',1);

 % Initialize variables that depend on number of nodes.
 neq = ndof*nnp;          % number of equations
 x = zeros(1,nnp);        % array of nodal x coordinates
 y = zeros(1,nnp);        % array of nodal y coordinates

 for i = 1:nnp
  % Read node id, but do not use it (it is assumed that
  % all nodes are numbered consecutively from one to the
  % total number of nodes).
  id = fscanf(fid,'%f',1);
  % Read nodal coordinates (disregard z coordinate).
  xyz = fscanf(fid,'%f',3);
  % armazena as coordenadas dos nós
  x(i) = xyz(1);
  y(i) = xyz(2);
 end
end

%-----------------------------------
function ReadNodeSupport(fid)
 include_gblrefs;

 % Store total number of nodes with essential B.C. (supports).
 nnebc = fscanf(fid,'%f',1);

 % Dimension arrays for support conditions.
 ebc_node = zeros(1,nnebc);     % array of nodes with essential B.C.
 ebc_flag = zeros(ndof,nnebc);  % array of nodal B.C. flags

 for i = 1:nnebc
  b = fscanf(fid,'%f',[1 7]);
   ebc_node(i) = b(1);          % get node id
   ebc_flag(:,i) = b(2:3);      % use only Dx and Dy constraints
 end
end

%-----------------------------------
function ReadMaterial(fid)
 include_gblrefs;

 % Store total number of materials.
 nmat = fscanf(fid,'%f',1);

 % Dimension vectors with material parameters.
 Elasticity = zeros(1,nmat); % array of material elasticity coeff values
 Poisson = zeros(1,nmat);    % array of material poisso ratio values
end

%-----------------------------------
function ReadMaterialIsotropic(fid)
 include_gblrefs;

 % Read number of materials in current file block.
 nmatb = fscanf(fid,'%f',1);

 for i = 1:nmatb
  % Read material id.
  id = fscanf(fid,'%f',1);
  % Read material properties parameters.
  mprops = fscanf(fid,'%f',2);
  % Store elasticity modulus and possion ratio of this material.
  if( id <= nmat )
   Elasticity(id) = mprops(1);
   Poisson(id) = mprops(2);
  end
 end
end

%-----------------------------------
function ReadThickness(fid)
 include_gblrefs;

 % Read number of specified thichnesses.
 nthk = fscanf(fid,'%f',1);
 Thickness = ones(nthk);

 for i = 1:nthk
  % Read thickness id (do not use it).
  id = fscanf(fid,'%f',1);
  % Read thickness value.
  Thickness(i) = fscanf(fid,'%f',1);
 end
end

%----------------------------------
function ReadIntegrationOrder(fid)
 global nintorder;                    % number of integration orders
 global intorder_data;                % integration order data

% Read integration order data.
 nintorder = fscanf(fid,'%f',1);
 intorder_data = zeros(nintorder,7);
 for i = 1:nintorder
  intorder_data(i,:) = fscanf(fid,'%f',[1 7]);
 end
end

%----------------------------------
function ReadElement(fid)
 include_gblrefs;

 % Read total number of elements in model.
 nel = fscanf(fid,'%f',1);

 % Dimension vectors of element attributes.
 emat = zeros(1,nel);          % array of element material id
 ethk = zeros(1,nel);          % array of element thickness id
end

%-----------------------------------
function ReadElementT3(fid)
 % 
 %                         3 +
 %                           |\
 %                           | \
 %                           |  \
 %                           |   \ T3 = TRIA3
 %                           |    \
 %                         1 +-----+ 2
 %
 include_gblrefs;
 global nintorder;                    % number of integration orders
 global intorder_data;                % integration order data

 element_type = TRIA3;
 edge_type = LINE2;
 nen = 3;                      % number of element nodes
 nedgen = 2;                   % number of edge nodes
 gauss_type = TRIA_QUADRATURE; % type of Gauss quadrature: triangle

 % Read number of elements in a neutral file block.
 nelem = fscanf(fid,'%f',1);
 for i = 1:nelem

  % Read current element data.
  b = fscanf(fid,'%f',[1 7]);

  % Skip element id (it is assumed that the elements are numbered
  % in consecutive order).

  % Store element material id.
  emat(i) = b(2);

  if(b(3) == 0)
   % Assume a default unit thickness property in case thickness id
   % is null.
   ethk(i) = 1;
  else
   % Store element thickness id.
   ethk(i) = b(3);
  end

  % Read id of integration order data.
  order_id = b(4);

  % Store element nodal incidence data.
  IEN(:,i) = b(5:7);
 end

  % Store integration order for element stiffness matrix computation.
  % Skip integration order id and use only the first direction order.
  % The last integration order id read is assumed as the valid one.
  gauss_order = intorder_data(order_id,2);
  % Store integration order for gauss point stress evaluation.
  % Use only the first direction order.
  gauss_stress_order = intorder_data(order_id,5);

end

%-----------------------------------
function ReadElementQ4(fid)
 % 
 %                         4 +---------------+ 3
 %                           |               |
 %                           |               |
 %                           |               |
 %                           |   Q4 = QUAD4  |
 %                           |               |
 %                           |               |
 %                         1 +---------------+ 2
 %
 include_gblrefs;
 global nintorder;                    % number of integration orders
 global intorder_data;                % integration order data

 element_type = QUAD4;
 edge_type = LINE2;
 nen = 4;                      % number of element nodes
 nedgen = 2;                   % number of edge nodes
 gauss_type = QUAD_QUADRATURE; % type of Gauss quadrature: quadrilateral

 % Read number of elements in a neutral file block.
 nelem = fscanf(fid,'%f',1);
 for i = 1:nelem

  % Read current element data.
  b = fscanf(fid,'%f',[1 8]);

  % Skip element id (it is assumed that the elements are numbered
  % in consecutive order).

  % Store element material id.
  emat(i) = b(2);

  if(b(3) == 0)
   % Assume a default unit thickness property in case thickness id
   % is null.
   ethk(i) = 1;
  else
   % Store element thickness id.
   ethk(i) = b(3);
  end

  % Read id of integration order data.
  order_id = b(4);

  % Store element nodal incidence data.
  IEN(:,i) = b(5:8);
 end

  % Store integration order for element stiffness matrix computation.
  % Skip integration order id and use only the first direction order.
  % The last integration order id read is assumed as the valid one.
  gauss_order = intorder_data(order_id,2);
  % Store integration order for gauss point stress evaluation.
  % Use only the first direction order.
  gauss_stress_order = intorder_data(order_id,5);

end

%-----------------------------------
function ReadElementT6(fid)
 % 
 %  Nodal numbering:     Neutral file        Current program
 %
 %                         5 +                3 +
 %                           |\                 |\
 %                           | \                | \
 %                         6 +  + 4           6 +  + 5
 %                           |   \ T6           |   \ TRIA6
 %                           |    \             |    \
 %                         1 +--+--+ 3        1 +--+--+ 2
 %                              2                  4
 %
 include_gblrefs;
 global nintorder;                    % number of integration orders
 global intorder_data;                % integration order data

 element_type = TRIA6;
 edge_type = LINE3;
 nen = 6;                      % number of element nodes
 nedgen = 3;                   % number of edge nodes
 gauss_type = TRIA_QUADRATURE; % type of Gauss quadrature: triangle

 % Read number of elements in a neutral file block.
 nelem = fscanf(fid,'%f',1);
 for i = 1:nelem

  % Read current element data.
  b = fscanf(fid,'%f',[1 10]);

  % Skip element id (it is assumed that the elements are numbered
  % in consecutive order).

  % Store element material id.
  emat(i) = b(2);

  if(b(3) == 0)
   % Assume a default unit thickness property in case thickness id
   % is null.
   ethk(i) = 1;
  else
   % Store element thickness id.
   ethk(i) = b(3);
  end

  % Read id of integration order data.
  order_id = b(4);

  % Store element nodal incidence data.
  IEN(1,i) = b(5);
  IEN(2,i) = b(7);
  IEN(3,i) = b(9);
  IEN(4,i) = b(6);
  IEN(5,i) = b(8);
  IEN(6,i) = b(10);
 end

  % Store integration order for element stiffness matrix computation.
  % Skip integration order id and use only the first direction order.
  % The last integration order id read is assumed as the valid one.
  gauss_order = intorder_data(order_id,2);
  % Store integration order for gauss point stress evaluation.
  % Use only the first direction order.
  gauss_stress_order = intorder_data(order_id,5);

end

%-----------------------------------
function ReadElementQ8(fid)
 % 
 % Nodal numbering:     Neutral file            Current program
 %
 %                          6                           7
 %                7 +-------+-------+ 5       4 +-------+-------+ 3
 %                  |               |           |               |
 %                  |               |           |               |
 %                8 +       Q8      + 4       8 +     QUAD8     + 6
 %                  |               |           |               |
 %                  |               |           |               |
 %                1 +-------+-------+ 3       1 +-------+-------+ 2
 %                          2                           5
 %
 include_gblrefs;
 global nintorder;                    % number of integration orders
 global intorder_data;                % integration order data

 element_type = QUAD8;
 edge_type = LINE3;
 nen = 8;                      % number of element nodes
 nedgen = 3;                   % number of edge nodes
 gauss_type = QUAD_QUADRATURE; % type of Gauss quadrature: quadrilateral

 % Read number of elements in a neutral file block.
 nelem = fscanf(fid,'%f',1);
 for i = 1:nelem

  % Read current element data.
  b = fscanf(fid,'%f',[1 12]);

  % Skip element id (it is assumed that the elements are numbered
  % in consecutive order).

  % Store element material id.
  emat(i) = b(2);

  if(b(3) == 0)
   % Assume a default unit thickness property in case thickness id
   % is null.
   ethk(i) = 1;
  else
   % Store element thickness id.
   ethk(i) = b(3);
  end

  % Read id of integration order data.
  order_id = b(4);

  % Store element nodal incidence data.
  IEN(1,i) = b(5);
  IEN(2,i) = b(7);
  IEN(3,i) = b(9);
  IEN(4,i) = b(11);
  IEN(5,i) = b(6);
  IEN(6,i) = b(8);
  IEN(7,i) = b(10);
  IEN(8,i) = b(12);
 end

  % Store integration order for element stiffness matrix computation.
  % Skip integration order id and use only the first direction order.
  % The last integration order id read is assumed as the valid one.
  gauss_order = intorder_data(order_id,2);
  % Store integration order for gauss point stress evaluation.
  % Use only the first direction order.
  gauss_stress_order = intorder_data(order_id,5);

end


%-----------------------------------
function ReadLoadCaseNodalForces(fid)
 include_gblrefs;

 % Read total number of point (node) loads.
 n_pointload = fscanf(fid,'%f',1);

 % Dimension arrays for nodal forces.
 pointload_node = zeros(1,n_pointload);      % array of nodes with point load
 pointload_value = zeros(ndof,n_pointload);  % array of point load values

 for i = 1:n_pointload
  b = fscanf(fid,'%f',[1 7]);
  % Store target node id and applied force component values (Fx and Fy).
  pointload_node(i) = b(1);
  pointload_value(:,i) = b(2:3);
 end
end

%-----------------------------------
function ReadLoadCaseNodalDisplacement(fid)
 include_gblrefs;

 % Read number of nodes with prescribed displacements.
 nnprescdispl = fscanf(fid,'%f',1);

 % Dimension arrays for prescribed displacements.
 prescdispl_node = zeros(1,nnprescdispl);     % array of nodes with presc.disp.
 prescdispl_value = zeros(ndof,nnprescdispl); % array of presc.disp. values

 for i = 1:nnprescdispl
  b = fscanf(fid,'%f',[1 7]);
  % Store target node id and prescribed displac. values (Dx and Dy).
  % Do not check for consistency with support condition nodes.
  prescdispl_node(i) = b(1);
  prescdispl_value(:,i) = b(2:3);
 end
end

%-----------------------------------
function ReadLoadCaseLineForceUniform(fid)
 include_gblrefs;

 % Read number of element sides with applied uniform distrib. loads.
 n_edgeload = fscanf(fid,'%f',1);   % number of edge uniform load
 for i=1:n_edgeload
  b = fscanf(fid,'%f',[1 7]);
  
  % Store id of element with applied edge uniform load.
  edgeload_side(1,i) = b(1);
  % Storte id's of nodes that define element edge.
  edgeload_side(2:3,i) = b(2:3);
  % Disregard information related to coordinate system: b(4).

  % Store uniform distributed edge load value.
  % (disregard z component)
  edgeload_value(:,i) = b(5:6);
 end
end

%-----------------------------------
function ReadLoadCaseDomainForceUniform(fid)
 include_gblrefs;

 % Read number of element with applied area uniform distrib. loads.
 n_areaload = fscanf(fid,'%f',1);   % number of area uniform load
 for i=1:n_areaload
  b = fscanf(fid,'%f',[1 4]);
  
  % Store id of element with applied area uniform load.
  areaload_elem(1,i) = b(1);

  % Store uniform distributed area load value.
  % (disregard z component)
  areaload_value(:,i) = b(5:6);
 end
end
