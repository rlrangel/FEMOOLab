function preReadNF(fid,model)
% This function reads a neutral-format file with input data for a finite-element model.
% The neutral format is described in www.tecgraf.puc-rio.br/neutralfile.
% It is assumed that all nodes and elements are numbered consecutively from
% one to the total number of nodes and elements.
% It is assumed that there is only one type of element and only one
% type of gauss integration order in each finite element model.
% It is assumed that there is a single load case.

flag = 1;
while flag
    % Get file line
    tline = fgetl(fid);
    % Get rid of blank spaces
    string = deblank(tline);
    
    % Check for a neutral-format tag string
    switch string
        case '%HEADER.ANALYSIS'
            ReadHeaderAnalysis(fid,model);
        case '%NODE.COORD'
            ReadNodeCoord(fid,model);
        case '%NODE.SUPPORT'
            ReadNodeSupport(fid,model);
        case '%MATERIAL'
            ReadMaterial(fid,model);
        case '%MATERIAL.ISOTROPIC'
            ReadMaterialIsotropic(fid,model);
        case '%THICKNESS'
            thickness = ReadThickness(fid);
        case '%INTEGRATION.ORDER'
            intorder_data = ReadIntegrationOrder(fid);
        case '%ELEMENT'
            ReadElement(fid,model);
        case '%ELEMENT.T3'
            ReadElementT3(fid,model,intorder_data,thickness);
        case '%ELEMENT.Q4'
            ReadElementQ4(fid,model,intorder_data,thickness);
        case '%ELEMENT.T6'
            ReadElementT6(fid,model,intorder_data,thickness);
        case '%ELEMENT.Q8'
            ReadElementQ8(fid,model,intorder_data,thickness);
        case '%LOAD.CASE.NODAL.FORCES'
            ReadLoadCaseNodalForces(fid,model);
        case '%LOAD.CASE.NODAL.DISPLACEMENT'
            ReadLoadCaseNodalDisplacement(fid,model);
        case '%LOAD.CASE.LINE.FORCE.UNIFORM'
            ReadLoadCaseLineForceUniform(fid,model);
        case '%LOAD.CASE.DOMAIN.FORCE.UNIFORM'
            ReadLoadCaseDomainForceUniform(fid,model);
    end
    
    if strcmp(string,'%END')
        flag = 0;
    end
end

fclose(fid);
end

%--------------------------------------------------------------------------
function ReadHeaderAnalysis(fid,model)
    % Get file line
    tline = fgetl(fid);
    % Get rid of blank spaces
    string = deblank(tline);

    % Get target analysis type
    switch string
        case '''plane_stress'''
            model.anm = Anm_PlaneStress();
        case '''plane_strain'''
            model.anm = Anm_PlaneStrain();
        case '''axisymmetric'''
            model.anm = Anm_Axisymmetric();
    end
end

%--------------------------------------------------------------------------
function ReadNodeCoord(fid,model)
    % Read total number of nodes in model
    nnp = fscanf(fid,'%f',1);
    model.nnp = nnp;
    model.neq = nnp * model.anm.ndof;

    % Initialize vector of objects of the Node class
    nodes(1,nnp) = Node();
    model.nodes = nodes;

    for n = 1:nnp
        b = fscanf(fid,'%f',4);
        % b(1) = node id
        % b(2) = X coordinate
        % b(3) = Y coordinate
        % b(4) = Z coordinate
        model.nodes(n).id = b(1);
        model.nodes(n).coords = [b(2) b(3) b(4)];
    end
end

%--------------------------------------------------------------------------
function ReadNodeSupport(fid,model)
    % Read total number of nodes with essential boundary conditions (supports)
    nnebc = fscanf(fid,'%f',1);
    
    for n = 1:nnebc
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Dx
        % b(3) = Dy
        % b(4) = Dz
        % b(5) = Rx
        % b(6) = Ry
        % b(7) = Rz
        model.nodes(b(1)).ebc = b(2:7);
    end
end

%--------------------------------------------------------------------------
function ReadMaterial(fid,model)
    % Read total number of materials in model
    nmat = fscanf(fid,'%f',1);
    model.nmat = nmat;
    
    % Initialize vector of objects of the Material class
    materials(1,nmat) = Material();
    model.materials = materials;
end

%--------------------------------------------------------------------------
function ReadMaterialIsotropic(fid,model)
    % Read number of materials in current file block
    nmatb = fscanf(fid,'%f',1);

    for m = 1:nmatb
        b = fscanf(fid,'%f',[1 3]);
        % b(1) = material id
        % b(2) = elasticity modulus
        % b(3) = Poisson ratio
        model.materials(m).id = b(1);
        if( b(1) <= model.nmat )
            model.materials(m).elasticity = b(2);
            model.materials(m).poisson = b(3);
        end
    end
end

%--------------------------------------------------------------------------
function thickness = ReadThickness(fid)
    % Read number of specified thicknesses
    nthk = fscanf(fid,'%f',1);
    thickness = ones(nthk,1);

    for t = 1:nthk
        b = fscanf(fid,'%f',[1 2]);
        % b(1) = thickness id
        % b(2) = thickness value
        thickness(t) = b(2);
    end
end

%--------------------------------------------------------------------------
function intorder_data = ReadIntegrationOrder(fid)
    % Read number of integration order data
    nintorder = fscanf(fid,'%f',1);
    intorder_data = zeros(nintorder,7);

    for i = 1:nintorder
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = ?
        % b(2) = ?
        % b(3) = ?
        % b(4) = ?
        % b(5) = ?
        % b(6) = ?
        % b(7) = ?
        intorder_data(i,:) = b;
    end
end

%--------------------------------------------------------------------------
function ReadElement(fid,model)
    % Read total number of elements in model
    model.nel = fscanf(fid,'%f',1);
end

%--------------------------------------------------------------------------
function ReadElementT3(fid,model,intorder_data,thickness)
%
%                         3 +
%                           |\
%                           | \
%                           |  \
%                           |   \ T3 = TRIA3
%                           |    \
%                         1 +-----+ 2
%
    % Read number of elements in a neutral file block
    nelem = fscanf(fid,'%f',1);
    
    % Initialize vector of objects of the Elem_T3 sub-class
    elems(1,nelem) = Elem_T3();
    model.elems = elems;

    for e = 1:nelem
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = element id
        % b(2) = element material id
        % b(3) = element thickness id
        % b(4) = id of integration order data
        % b(5:7) = element nodal incidence data
        model.elems(e).id = b(1);
        model.elems(e).material = model.materials(b(2));
        if b(3) == 0
            model.elems(e).thk = 1;
        else
            model.elems(e).thk = thickness(b(3));
        end
        model.elems(e).gauss_order = intorder_data(b(4),2);
        model.elems(e).nodes = [ model.nodes(b(5)) model.nodes(b(6)) model.nodes(b(7)) ];
        
        % Set element analysis model
        model.elems(e).anm = model.anm;
    end

    % Store integration order for gauss point stress evaluation.
    % Use only the first direction order.
    %gauss_stress_order = intorder_data(order_id,5);
end

%--------------------------------------------------------------------------
function ReadElementQ4(fid,model,intorder_data,thickness)
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
    % Read number of elements in a neutral file block
    nelem = fscanf(fid,'%f',1);
    
    % Initialize vector of objects of the Elem_Q4 sub-class
    elems(1,nelem) = Elem_Q4();
    model.elems = elems;

    for e = 1:nelem
        b = fscanf(fid,'%f',[1 8]);
        % b(1) = element id
        % b(2) = element material id
        % b(3) = element thickness id
        % b(4) = id of integration order data
        % b(5:8) = element nodal incidence data
        model.elems(e).id = b(1);
        model.elems(e).material = model.materials(b(2));
        model.elems(e).thk = thickness(b(3));
        if b(3) == 0
            model.elems(e).thk = 1;
        else
            model.elems(e).thk = thickness(b(3));
        end
        model.elems(e).gauss_order = intorder_data(b(4),2);
        model.elems(e).nodes = [ model.nodes(b(5)) model.nodes(b(6))...
                                 model.nodes(b(7)) model.nodes(b(8)) ];
        
        % Set element analysis model
        model.elems(e).anm = model.anm;
    end

    % Store integration order for gauss point stress evaluation.
    % Use only the first direction order.
    %gauss_stress_order = intorder_data(order_id,5);
end

%--------------------------------------------------------------------------
function ReadElementT6(fid,model,intorder_data,thickness)
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
    % Read number of elements in a neutral file block
    nelem = fscanf(fid,'%f',1);
    
    % Initialize vector of objects of the Elem_T6 sub-class
    elems(1,nelem) = Elem_T6();
    model.elems = elems;

    for e = 1:nelem
        b = fscanf(fid,'%f',[1 10]);
        % b(1) = element id
        % b(2) = element material id
        % b(3) = element thickness id
        % b(4) = id of integration order data
        % b(5:10) = element nodal incidence data
        model.elems(e).id = b(1);
        model.elems(e).material = model.materials(b(2));
        model.elems(e).thk = thickness(b(3));
        if b(3) == 0
            model.elems(e).thk = 1;
        else
            model.elems(e).thk = thickness(b(3));
        end
        model.elems(e).gauss_order = intorder_data(b(4),2);
        model.elems(e).nodes = [ model.nodes(b(5)) model.nodes(b(7)) model.nodes(b(9))...
                                 model.nodes(b(6)) model.nodes(b(8)) model.nodes(b(10)) ];
        
        % Set element analysis model
        model.elems(e).anm = model.anm;
    end

    % Store integration order for gauss point stress evaluation.
    % Use only the first direction order.
    %gauss_stress_order = intorder_data(order_id,5);
end

%--------------------------------------------------------------------------
function ReadElementQ8(fid,model,intorder_data,thickness)
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
    % Read number of elements in a neutral file block
    nelem = fscanf(fid,'%f',1);
    
    % Initialize vector of objects of the Elem_Q8 sub-class
    elems(1,nelem) = Elem_Q8();
    model.elems = elems;

    for e = 1:nelem
        b = fscanf(fid,'%f',[1 12]);
        % b(1) = element id
        % b(2) = element material id
        % b(3) = element thickness id
        % b(4) = id of integration order data
        % b(5:12) = element nodal incidence data
        model.elems(e).id = b(1);
        model.elems(e).material = model.materials(b(2));
        model.elems(e).thk = thickness(b(3));
        if b(3) == 0
            model.elems(e).thk = 1;
        else
            model.elems(e).thk = thickness(b(3));
        end
        model.elems(e).gauss_order = intorder_data(b(4),2);
        model.elems(e).nodes = [ model.nodes(b(5))  model.nodes(b(7)) ...
                                 model.nodes(b(9))  model.nodes(b(11))...
                                 model.nodes(b(6))  model.nodes(b(8)) ...
                                 model.nodes(b(10)) model.nodes(b(12)) ];
        
        % Set element analysis model
        model.elems(e).anm = model.anm;
    end

    % Store integration order for gauss point stress evaluation.
    % Use only the first direction order.
    %gauss_stress_order = intorder_data(order_id,5);
end

%--------------------------------------------------------------------------
function ReadLoadCaseNodalForces(fid,model)
    % Read total number of point (node) loads
    n_pointload = fscanf(fid,'%f',1);

    for i = 1:n_pointload
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Fx
        % b(3) = Fy
        % b(4) = Fz
        % b(5) = Mx
        % b(6) = My
        % b(7) = Mz
        model.nodes(b(1)).nodalLoad = b(2:7);
    end
end

%--------------------------------------------------------------------------
function ReadLoadCaseNodalDisplacement(fid,model)
    % Read number of nodes with prescribed displacements
    nnprescdispl = fscanf(fid,'%f',1);

    for i = 1:nnprescdispl
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Dx
        % b(3) = Dy
        % b(4) = Dz
        % b(5) = Rx
        % b(6) = Ry
        % b(7) = Rz
        model.nodes(b(1)).prescDispl = b(2:7);
    end
end

%--------------------------------------------------------------------------
function ReadLoadCaseLineForceUniform(fid,model)
    % Read number of element sides with applied uniform distrib. loads
    n_edgeload = fscanf(fid,'%f',1);
    
    a = zeros(n_edgeload,1);
    for i = 1:n_edgeload
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = element id
        % b(2) = initial node id
        % b(3) = final node id
        % b(4) = coordinate system [0=global, 1=local]
        % b(5) = qx
        % b(6) = qy
        % b(7) = qz
        a(i) = b(1);
        l = sum(a(:)==b(1));
        model.elems(b(1)).edgeLoad(l,:) = [ b(2) b(3) b(5) b(6) b(7) ];
    end
end

%--------------------------------------------------------------------------
function ReadLoadCaseDomainForceUniform(fid,model)
    % Read number of element with applied area uniform distrib. loads
    n_areaload = fscanf(fid,'%f',1);

    for i=1:n_areaload
        b = fscanf(fid,'%f',[1 8]);
        % b(1) = element id
        % b(2) = node i id
        % b(3) = node j id
        % b(4) = node k id
        % b(5) = coordinate system [0=global, 1=local]
        % b(6) = px
        % b(7) = py
        % b(8) = pz
        model.elems(b(1)).areaLoad = b(6:8);
    end
end
