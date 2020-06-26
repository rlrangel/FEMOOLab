%% Read Class
%
% This class defines a reader object in the StAnOOP program.
% A reader object is responsible for reading the input file with model/analysis
% information and store it in the program data structure.
%
classdef Read < handle
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function read = Read()
            return;
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Main function to read input file tags. If file is read successfully,
        % status is 1, otherwise it is 0.
        function status = inputFile(read,fin,sim)
            status = 1;
            
            % Create Model, Analysis and Result objects
            sim.mdl = Model();
            sim.anl = Anl_LinearElastic();
            sim.res = Result();
            
            while status == 1 && ~feof(fin)
                % Get file line without blank spaces
                string = deblank(fgetl(fin));
                
                % Look for tag strings
                switch string
                    case '%MODEL.ANALYSIS_MODEL'
                        status = read.analysisModel(fin,sim.mdl);
                    case '%MODEL.NODE_COORD'
                        status = read.nodeCoord(fin,sim.mdl);
                    case '%MODEL.NODE_SUPPORT'
                        status = read.nodeSupport(fin,sim.mdl);
                    case '%MODEL.NODE_PRESC_DISPL'
                        status = read.nodePrescDispl(fin,sim.mdl);
                    case '%MODEL.MATERIAL_PROPERTY'
                        status = read.materialProperty(fin,sim.mdl);
                    case '%MODEL.ELEMENT_TRI3'
                        status = read.elementTri3(fin,sim.mdl);
                    case '%MODEL.ELEMENT_QUAD4'
                        status = read.elementQuad4(fin,sim.mdl);
                    case '%MODEL.ELEMENT_TRI6'
                        status = read.elementTri6(fin,sim.mdl);
                    case '%MODEL.ELEMENT_QUAD8'
                        status = read.elementQuad8(fin,sim.mdl);
                    case '%MODEL.LOAD_POINT'
                        status = read.loadPoint(fin,sim.mdl);
                    case '%MODEL.LOAD_LINE_UNIFORM'
                        status = read.loadLineUnif(fin,sim.mdl);
                    case '%MODEL.LOAD_AREA_UNIFORM'
                        status = read.loadAreaUnif(fin,sim.mdl);
                end
                
                if (strcmp(string,'%END'))
                    break;
                end
            end
            if (status == 1)
                status = read.checkInput(sim.mdl);
            end
            fclose(fin);
        end
    
    %% Methods for reading TAGS
        %------------------------------------------------------------------
        function status = analysisModel(~,fin,mdl)
            status = 1;
            string = deblank(fgetl(fin));
            
            if (strcmp(string,'''PLANE_STRESS'''))
                mdl.anm = Anm_PlaneStress();
            elseif (strcmp(string,'''PLANE_STRAIN'''))
                mdl.anm = Anm_PlaneStrain();
            elseif (strcmp(string,'''AXISYMMETRIC'''))
                mdl.anm = Anm_Axisymmetric();
            else
                fprintf('Invalid input data: analysis model!\n');
                status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function status = nodeCoord(read,fin,mdl)
            status = 1;
            if (isempty(mdl.anm))
                fprintf('Analysis model must be provided before node coordinates!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of nodes'))
                status = 0;
                return;
            end
            
            % Create vector of Node objects and compute total number of equations
            mdl.nnp = n;
            nodes(n,1) = Node();
            mdl.nodes = nodes;
            mdl.neq = n * mdl.anm.ndof;
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,n,'node ID for coordinate specification'))
                    status = 0;
                    return;
                end
                
                % Coordinates (X,Y,Z)
                [coord,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid coordinates of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).id = id;
                mdl.nodes(id).coord = coord';
            end
        end
        
        %------------------------------------------------------------------
        function status = nodeSupport(read,fin,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before node supports!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with supports
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,mdl.nnp,'number of nodes with support conditions'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nnp,'node ID for support condition specification'))
                    status = 0;
                    return;
                end
                
                % Support condition flags
                [ebc,count] = fscanf(fin,'%d',3);
                if (count~= 3 || ~all(ismember(ebc,0:1)))
                    fprintf('Invalid support conditions of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).ebc = ebc';
            end
        end
        
        %------------------------------------------------------------------
        function status = nodePrescDispl(read,fin,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before prescribed displacements!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with prescribed displacements
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,mdl.nnp,'number of nodes with prescribed displacements'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nnp,'node ID for prescribed displacement specification'))
                    status = 0;
                    return;
                end
                
                % Prescribed displacements values
                [disp,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid prescribed displacements of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).prescDispl = disp';
            end
        end
        
        %------------------------------------------------------------------
        function status = materialProperty(read,fin,mdl)
            status = 1;
            
            % Total number of materials
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of materials'))
                status = 0;
                return;
            end
            
            % Create vector of Material objects
            mdl.nmat = n;
            materials(n,1) = Material();
            mdl.materials = materials;
            
            for i = 1:n
                % Material ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,n,'material ID'))
                    status = 0;
                    return;
                end
                
                % Material properties: E,v
                [prop,count] = fscanf(fin,'%f',2);
                if (~read.chkMatProp(prop,count,id))
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.materials(id).id = id;
                mdl.materials(id).E  = prop(1);
                mdl.materials(id).v  = prop(2);
            end
        end
        
        %------------------------------------------------------------------
        function status = elementTri3(read,fin,mdl)
            %                         3 +
            %                           |\
            %                           | \
            %                           |  \
            %                           |   \ TRI3
            %                           |    \
            %                         1 +-----+ 2
            status = 1;
            if (isempty(mdl.nodes) || isempty(mdl.materials))
                fprintf('Node coordinates and materials must be provided before elements!\n');
                status = 0;
                return;
            end
            
            % Number of elements T3
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of elements'))
                status = 0;
                return;
            end
            
            % Create vector of Element Objects
            mdl.nel = n;
            elems(n,1) = Element_Tri3();
            mdl.elems = elems;
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nel,'element ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Thickness, Material ID, Gauss quadrature orders, Connectivity
                [p,count] = fscanf(fin,'%f',7);
                if (~read.chkElemProp(mdl,id,p,count,7))
                    status = 0;
                    return;
                end
                
                mdl.elems(id) = Element_Tri3(...
                                id,mdl.anm,p(1),mdl.materials(p(2)),p(3),p(4),...
                                [mdl.nodes(p(5));mdl.nodes(p(6));mdl.nodes(p(7))]);
            end
        end
        
        %------------------------------------------------------------------
        function status = elementQuad4(read,fin,mdl)
            %                     4 +---------------+ 3
            %                       |               |
            %                       |               |
            %                       |     QUAD      |
            %                       |               |
            %                       |               |
            %                     1 +---------------+ 2
            status = 1;
            if (isempty(mdl.nodes) || isempty(mdl.materials))
                fprintf('Node coordinates and materials must be provided before elements!\n');
                status = 0;
                return;
            end
            
            % Number of elements Q4
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of elements'))
                status = 0;
                return;
            end
            
            % Create vector of Element Objects
            mdl.nel = n;
            elems(n,1) = Element_Quad4();
            mdl.elems = elems;
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nel,'element ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Thickness, Material ID, Gauss quadrature orders, Connectivity
                [p,count] = fscanf(fin,'%f',8);
                if (~read.chkElemProp(mdl,id,p,count,8))
                    status = 0;
                    return;
                end
                
                mdl.elems(id) = Element_Quad4(...
                                id,mdl.anm,p(1),mdl.materials(p(2)),p(3),p(4),...
                                [mdl.nodes(p(5));mdl.nodes(p(6));...
                                 mdl.nodes(p(7));mdl.nodes(p(8))]);
            end
        end
        
        %------------------------------------------------------------------
        function status = elementTri6(read,fin,mdl)
            %                         3 +
            %                           |\
            %                           | \
            %                         6 +  + 5
            %                           |   \ TRI6
            %                           |    \
            %                         1 +--+--+ 2
            %                              4
            status = 1;
            if (isempty(mdl.nodes) || isempty(mdl.materials))
                fprintf('Node coordinates and materials must be provided before elements!\n');
                status = 0;
                return;
            end
            
            % Number of elements T6
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of elements'))
                status = 0;
                return;
            end
            
            % Create vector of Element Objects
            mdl.nel = n;
            elems(n,1) = Element_Tri6();
            mdl.elems = elems;
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nel,'element ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Thickness, Material ID, Gauss quadrature orders, Connectivity
                [p,count] = fscanf(fin,'%f',10);
                if (~read.chkElemProp(mdl,id,p,count,10))
                    status = 0;
                    return;
                end
                
                mdl.elems(id) = Element_Tri6(...
                                id,mdl.anm,p(1),mdl.materials(p(2)),p(3),p(4),...
                                [mdl.nodes(p(5));mdl.nodes(p(6));...
                                 mdl.nodes(p(7));mdl.nodes(p(8));...
                                 mdl.nodes(p(9));mdl.nodes(p(10))]);
            end
        end
        
        %------------------------------------------------------------------
        function status = elementQuad8(read,fin,mdl)
            %                              7
            %                    4 +-------+-------+ 3
            %                      |               |
            %                      |               |
            %                    8 +     QUAD8     + 6
            %                      |               |
            %                      |               |
            %                    1 +-------+-------+ 2
            %                              5
            status = 1;
            if (isempty(mdl.nodes) || isempty(mdl.materials))
                fprintf('Node coordinates and materials must be provided before elements!\n');
                status = 0;
                return;
            end
            
            % Number of elements Q8
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of elements'))
                status = 0;
                return;
            end
            
            % Create vector of Element Objects
            mdl.nel = n;
            elems(n,1) = Element_Quad8();
            mdl.elems = elems;
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nel,'element ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Thickness, Material ID, Gauss quadrature orders, Connectivity
                [p,count] = fscanf(fin,'%f',12);
                if (~read.chkElemProp(mdl,id,p,count,12))
                    status = 0;
                    return;
                end
                
                mdl.elems(id) = Element_Quad8(...
                                id,mdl.anm,p(1),mdl.materials(p(2)),p(3),p(4),...
                                [mdl.nodes(p(5)); mdl.nodes(p(6));...
                                 mdl.nodes(p(7)); mdl.nodes(p(8));...
                                 mdl.nodes(p(9)); mdl.nodes(p(10));...
                                 mdl.nodes(p(11));mdl.nodes(p(12))]);
            end
        end
        
        %--------------------------------------------------------------------------
        function status = loadPoint(read,fin,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before point loads!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with point load
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,mdl.nnp,'number of nodes with load'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nnp,'node ID for point load specification'))
                    status = 0;
                    return;
                end
                
                % Point load values
                [load,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid point load of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).prescDispl = load';
            end
        end
        
        %--------------------------------------------------------------------------
        function status = loadLineUnif(read,fin,mdl)
            status = 1;
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before line loads!\n');
                status = 0;
                return;
            end
            
            % Total number of edges with line load
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of edges with line load'))
                status = 0;
                return;
            end
            
            a = zeros(n,1);
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nel,'element ID for line load specification'))
                    status = 0;
                    return;
                end
                
                % Area load information (node1,node2,qx,qy,qz)
                [load,count] = fscanf(fin,'%f',5);
                if (count ~= 5)
                    fprintf('Invalid line load specification of element %d\n',id);
                    status = 0;
                    return;
                end
                % CHECK IF NODES ARE ON THE SAME EDGE !!!
                
                % Store data
                a(i) = id;
                j = sum(a(:)==id);
                mdl.elems(id).lineLoad(j,:) = load';
            end
        end
        
        %--------------------------------------------------------------------------
        function status = loadAreaUnif(read,fin,mdl)
            status = 1;
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before area loads!\n');
                status = 0;
                return;
            end
            
            % Total number of edges with line load
            n = fscanf(fin,'%d',1);
            if (~read.chkInt(n,inf,'number of faces with area load'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                if (~read.chkInt(id,mdl.nel,'element ID for area load specification'))
                    status = 0;
                    return;
                end
                
                % Area load information (qx,qy,qz)
                [load,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid area load specification of element %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.elems(id).areaLoad = load';
            end
        end
    
    %% Static methods for checking input data
        %------------------------------------------------------------------
        function status = checkInput(~,mdl)
            status = 0;
            if isempty(mdl.anm)
                fprintf('Analysis model type not provided!\n');
                return;
            end
            if isempty(mdl.nodes)
                fprintf('Nodes not provided!\n');
                return;
            end
            if isempty(mdl.elems)
                fprintf('Elements not provided!\n');
                return;
            end
            if isempty(mdl.materials)
                fprintf('Materials not provided!\n');
                return;
            end
            status = 1;
        end
        
        %------------------------------------------------------------------
        function status = chkInt(~,val,max,string)
            status = 1;
            if (~isnumeric(val)         ||...
                floor(val) ~= ceil(val) ||...
                val <= 0                ||...
                val > max)
                string = strcat('Invalid: ',string,'!\n');
                fprintf(string);
                if (max == inf)
                    fprintf('It must be a positive integer!\n');
                else
                    fprintf('It must be a positive integer less/equal to %d!\n',max);
                end
                status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function status = chkMatProp(~,prop,count,id)
            status = 0;
            if (count ~= 2)
                fprintf('Invalid properties for material %d\n',id);
            elseif (prop(1) <= 0)
                fprintf('Invalid properties for material %d\n',id);
                fprintf('Elasticity modulus must be positive!');
            elseif (prop(2) <= -1 || prop(2) > 0.5)
                fprintf('Invalid properties for material %d\n',id);
                fprintf('Poisson ratio must be between -1.0 and +0.5!');
            else
                status = 1;
            end
        end
        
        %------------------------------------------------------------------
        function status = chkElemProp(read,mdl,id,prop,count,num)
            status = 0;
            nodesIDs = prop(5:end);
            if (count ~= num)
                fprintf('Invalid properties of element %d\n',id);
            elseif (prop(1) <= 0)
                fprintf('Invalid properties of element %d\n',id);
                fprintf('Thickness must be positive!\n');
            elseif (~read.chkInt(prop(2),mdl.nmat,'material ID for element definition'))
            elseif (~read.chkInt(prop(3),4,'Gauss quadrature order for element definition'))
            elseif (~read.chkInt(prop(4),4,'Gauss quadrature order for element definition'))
            elseif (min(nodesIDs) <= 0 || max(nodesIDs) > mdl.nnp)
                fprintf('Invalid properties of element %d\n',id);
                fprintf('Node IDs must be a positive integer less/equal to the total number of nodes!\n');
            elseif (length(nodesIDs) ~= length(unique(nodesIDs)))
                fprintf('Invalid properties of element %d\n',id);
                fprintf('Repeated node ID!\n');
            else
                status = 1;
            end
        end
    end
end