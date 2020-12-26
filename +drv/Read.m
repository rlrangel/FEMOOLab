%% Read Class
%
%% Description
%
% This class defines a reader object in the FEMOOLab program.
% A reader object is responsible for reading a
% <https://web.tecgraf.puc-rio.br/neutralfile neutral file> with model and
% analysis information and store it in the simulation data structure.
%
%% Class definition
%
classdef Read < handle
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Read(fid,sim,opt)
            if (nargin == 3)
                this.execute(fid,sim,opt);
            end
        end
    end
    
    %% Public methods
    methods
        %% Main function
        %------------------------------------------------------------------
        function status = execute(this,fid,sim,opt)
            status = 1;
            
            % Set result plotting options (they will be read from file in the future)
            this.setResultPlottingOpt(sim.mdl,opt);
            
            % Create Gauss quadrature for triangular and quadrilateral shapes
            gauss_tria = fem.Gauss_Tria();
            gauss_quad = fem.Gauss_Quad();
            
            while status == 1 && ~feof(fid)
                % Get file line without blank spaces
                string = deblank(fgetl(fid));
                
                % Look for tag strings
                switch string
                    case '%HEADER.ANALYSIS'
                        status = this.analysisModel(fid,sim);
                    case '%NODE.COORD'
                        status = this.nodeCoord(fid,sim.mdl);
                    case '%NODE.SUPPORT'
                        status = this.nodeSupport(fid,sim.mdl);
                    case '%MATERIAL'
                        status = this.materialTotal(fid,sim.mdl);
                    case '%MATERIAL.ISOTROPIC'
                        status = this.materialIsotropic(fid,sim.mdl);
                    case '%MATERIAL.PROPERTY.DENSITY'
                        status = this.materialDensity(fid,sim.mdl);
                    case '%MATERIAL.PROPERTY.THERMAL'
                        status = this.materialThermal(fid,sim.mdl);
                    case '%THICKNESS'
                        [status,thickness] = this.Thickness(fid);
                    case '%INTEGRATION.ORDER'
                        [status,intgrorder] = this.IntgrOrder(fid);
                    case '%ELEMENT'
                        status = this.elementTotal(fid,sim.mdl);
                    case '%ELEMENT.T3'
                        status = this.elementTria3(fid,sim.mdl,gauss_tria,thickness,intgrorder);
                    case '%ELEMENT.Q4'
                        status = this.elementQuad4(fid,sim.mdl,gauss_quad,thickness,intgrorder);
                    case '%ELEMENT.T6'
                        status = this.elementTria6(fid,sim.mdl,gauss_tria,thickness,intgrorder);
                    case '%ELEMENT.Q8'
                        status = this.elementQuad8(fid,sim.mdl,gauss_quad,thickness,intgrorder);
                    case '%LOAD.CASE.NODAL.FORCES'
                        status = this.loadPoint(fid,sim.mdl);
                    case '%LOAD.CASE.NODAL.DISPLACEMENT'
                        status = this.nodePrescDispl(fid,sim.mdl);
                    case '%LOAD.CASE.NODAL.TEMPERATURE'
                        status = this.nodePrescTemp(fid,sim.mdl);
                    case '%LOAD.CASE.LINE.FORCE.UNIFORM'
                        status = this.loadLineUnif(fid,sim.mdl);
                    case '%LOAD.CASE.DOMAIN.FORCE.UNIFORM'
                        status = this.loadDomainUnif(fid,sim.mdl);
                    case '%LOAD.CASE.LINE.HEAT.FLUX.UNIFORM'
                        status = this.fluxLineUnif(fid,sim.mdl);
                    case '%LOAD.CASE.AREA.HEAT.FLUX.UNIFORM' % This is actually DOMAIN heat flux!
                        status = this.fluxDomainUnif(fid,sim.mdl);
                end
                if (strcmp(string,'%END'))
                    break;
                end
            end
            if (status == 1)
                status = this.checkInput(sim.mdl);
            end
            fclose(fid);
        end
        
        %% Method for reading plotting options (temporary)
        %------------------------------------------------------------------
        function setResultPlottingOpt(~,mdl,opt)
            mdl.res.eid    = opt.eid;
            mdl.res.nid    = opt.nid;
            mdl.res.scl    = opt.scl;
            mdl.res.dx     = opt.dx;
            mdl.res.dy     = opt.dy;
            mdl.res.dz     = opt.dz;
            mdl.res.rx     = opt.rx;
            mdl.res.ry     = opt.ry;
            mdl.res.rz     = opt.rz;
            mdl.res.temp   = opt.temp;
            mdl.res.smooth = opt.smooth;
            mdl.res.sxx    = opt.sxx;
            mdl.res.syy    = opt.syy;
            mdl.res.szz    = opt.szz;
            mdl.res.txy    = opt.txy;
            mdl.res.txz    = opt.txz;
            mdl.res.tyz    = opt.tyz;
            mdl.res.s1     = opt.s1;
            mdl.res.s2     = opt.s2;
            mdl.res.s3     = opt.s3;
            mdl.res.taumax = opt.taumax;
            mdl.res.mxx    = opt.mxx;
            mdl.res.myy    = opt.myy;
            mdl.res.mxy    = opt.mxy;
            mdl.res.qxz    = opt.qxz;
            mdl.res.qyz    = opt.qyz;
            mdl.res.m1     = opt.m1;
            mdl.res.m2     = opt.m2;
            mdl.res.tormax = opt.tormax;
            mdl.res.fxx    = opt.fxx;
            mdl.res.fyy    = opt.fyy;
            mdl.res.fzz    = opt.fzz;
        end
        
        %% Methods for reading neutral file TAGS
        %------------------------------------------------------------------
        function status = analysisModel(~,fid,sim)
            status = 1;
            string = deblank(fgetl(fid));
            
            if (strcmp(string,'''PLANE_STRESS''') || strcmp(string,'''plane_stress'''))
                sim.mdl.anm = fem.Anm_PlaneStress();
                sim.anl = fem.Anl_LinearStatic();
            elseif (strcmp(string,'''PLANE_STRAIN''') || strcmp(string,'''plane_strain'''))
                sim.mdl.anm = fem.Anm_PlaneStrain();
                sim.anl = fem.Anl_LinearStatic();
            elseif (strcmp(string,'''PLANE_CONDUCTION''') || strcmp(string,'''plane_conduction'''))
                sim.mdl.anm = fem.Anm_PlaneConduction();
                sim.anl = fem.Anl_LinearStatic();
            elseif (strcmp(string,'''AXISYMMETRIC''') || strcmp(string,'''axisymmetric'''))
                sim.mdl.anm = fem.Anm_Axisymmetric();
                sim.anl = fem.Anl_LinearStatic();
            else
                fprintf('Invalid input data: analysis model!\n');
                status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function status = nodeCoord(this,fid,mdl)
            status = 1;
            if (isempty(mdl.anm))
                fprintf('Analysis model must be provided before node coordinates!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of nodes'))
                status = 0;
                return;
            end
            
            % Create vector of Node objects
            mdl.nnp = n;
            nodes(n,1) = fem.Node();
            mdl.nodes = nodes;
            
            for i = 1:n
                % Node ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,n,'node ID for coordinate specification'))
                    status = 0;
                    return;
                end
                
                % Coordinates (X,Y,Z)
                [coord,count] = fscanf(fid,'%f',3);
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
        function status = nodeSupport(this,fid,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before node supports!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with supports
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nnp,'number of nodes with support conditions'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nnp,'node ID for support condition specification'))
                    status = 0;
                    return;
                end
                
                % Support condition flags
                [supp,count] = fscanf(fid,'%d',6);
                if (count ~= 6 || ~all(ismember(supp,0:1)))
                    fprintf('Invalid support conditions of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).fixDispl = logical(supp);
            end
        end
        
        %------------------------------------------------------------------
        function status = materialTotal(this,fid,mdl)
            status = 1;

            % Total number of materials
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'total number of materials'))
                status = 0;
                return;
            end
            
            % Create vector of Material objects
            mdl.nmat = n;
            materials(n,1) = fem.Material();
            mdl.materials = materials;
            
            % Set materials IDs
            for i = 1:n
                mdl.materials(i).id = i;
            end
        end
        
        %------------------------------------------------------------------
        function status = materialIsotropic(this,fid,mdl)
            status = 1;
            if (isempty(mdl.materials))
                fprintf('Total number of materials must be provided before isotropic properties!\n');
                status = 0;
                return;
            end
            
            % Number of material isotropic properties
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nmat,'number of material isotropic properties'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Material ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nmat,'material ID for isotropic properties'))
                    status = 0;
                    return;
                end
                
                % Material properties: E,v
                [prop,count] = fscanf(fid,'%f',2);
                if (count ~= 2)
                    fprintf('Invalid isotropic properties for material %d\n',id);
                    status = 0;
                    return;
                elseif (prop(1) <= 0)
                    fprintf('Invalid isotropic properties for material %d\n',id);
                    fprintf('Elasticity modulus must be positive!');
                    status = 0;
                    return;
                elseif (prop(2) <= -1 || prop(2) > 0.5)
                    fprintf('Invalid isotropic properties for material %d\n',id);
                    fprintf('Poisson ratio must be between -1.0 and +0.5!');
                    status = 0;
                    return;
                end
                
                % Store isotropic properties
                mdl.materials(id).E  = prop(1);
                mdl.materials(id).v  = prop(2);
            end
        end
        
        %------------------------------------------------------------------
        function status = materialDensity(this,fid,mdl)
           status = 1;
            if (isempty(mdl.materials))
                fprintf('Total number of materials must be provided before densities!\n');
                status = 0;
                return;
            end
            
            % Number of material densities
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nmat,'number of material densities'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Material ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nmat,'material ID for densities'))
                    status = 0;
                    return;
                end
                
                % Material densities
                [rho,count] = fscanf(fid,'%f',1);
                if (count ~= 1)
                    fprintf('Invalid density for material %d\n',id);
                    status = 0;
                    return;
                elseif (rho <= 0)
                    fprintf('Invalid properties for material %d\n',id);
                    fprintf('Density must be positive!');
                    status = 0;
                    return;
                end
                
                % Store density
                mdl.materials(id).rho = rho;
            end
        end
        
        %------------------------------------------------------------------
        function status = materialThermal(this,fid,mdl)
            status = 1;
            if (isempty(mdl.materials))
                fprintf('Total number of materials must be provided before thermal properties!\n');
                status = 0;
                return;
            end
            
            % Number of material thermal properties
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nmat,'number of material thermal properties'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Material ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nmat,'material ID for thermal properties'))
                    status = 0;
                    return;
                end
                
                % Material properties: cp,k
                [prop,count] = fscanf(fid,'%f',2);
                if (count ~= 2)
                    fprintf('Invalid thermal properties for material %d\n',id);
                    status = 0;
                    return;
                elseif (prop(1) <= 0)
                    fprintf('Invalid thermal properties for material %d\n',id);
                    fprintf('Thermal conductivity must be positive!');
                    status = 0;
                    return;
                elseif (prop(2) <= 0)
                    fprintf('Invalid thermal properties for material %d\n',id);
                    fprintf('Specific heat capacity must be positive!');
                    status = 0;
                    return;
                end
                
                % Store thermal properties
                mdl.materials(id).k  = prop(1);
                mdl.materials(id).cp = prop(2);
            end
        end
        
        %------------------------------------------------------------------
        function [status,thickness] = Thickness(this,fid)
            status = 1;
            
            % Total number of thickness values
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of thicknesses'))
                status = 0;
                return;
            end
            
            % Create vector of thicknesses
            thickness = zeros(n,1);
            
            for i = 1:n
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,n,'ID of thickness'))
                    status = 0;
                    return;
                end
                
                [thk,count] = fscanf(fid,'%f',1);
                if (count ~= 1)
                    fprintf('Invalid thickness %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                thickness(id) = thk;
            end
        end
        
        %------------------------------------------------------------------
        function [status,intgrorder] = IntgrOrder(this,fid)
            status = 1;
            
            % Total number of integration orders
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of integration orders'))
                status = 0;
                return;
            end

            % Create vector of integration orders
            intgrorder = zeros(n,2);
            
            for i = 1:n
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,n,'ID of integration order'))
                    status = 0;
                    return;
                end
                
                % Stiffness and stress integration orders
                [order,count] = fscanf(fid,'%d',6);
                if (count ~= 6)
                    fprintf('Invalid integration order %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                intgrorder(id,1) = order(1);
                intgrorder(id,2) = order(4);
            end
        end
        
        %------------------------------------------------------------------
        function status = elementTotal(this,fid,mdl)
            status = 1;

            % Total number of elements
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of elements'))
                status = 0;
                return;
            end
            
            % Create vector of Element objects
            mdl.nel = n;
            elems(n,1) = fem.Element();
            mdl.elems = elems;
            
            % Set elements IDs
            for i = 1:n
                mdl.elems(i).id = i;
            end
        end
        
        %------------------------------------------------------------------
        function status = elementTria3(this,fid,mdl,gauss,thickness,intgrorder)
            %                         3 +
            %                           |\
            %                           | \
            %                           |  \
            %                           |   \ TRIA3
            %                           |    \
            %                         1 +-----+ 2
            status = 1;
            if (isempty(mdl.nodes) || isempty(mdl.materials))
                fprintf('Node coordinates and materials must be provided before elements!\n');
                status = 0;
                return;
            end
            
            % Number of elements T3
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nel,'number of elements T3'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element T3 ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Material ID, Thickness ID, Integration order ID, Connectivity
                [p,count] = fscanf(fid,'%d',6);
                if (~this.chkElemProp(mdl,id,p,thickness,intgrorder,count,6))
                    status = 0;
                    return;
                end
                
                % Create shape object
                shape = fem.Shape_Tria3([mdl.nodes(p(4));mdl.nodes(p(5));mdl.nodes(p(6))]);
                
                % Store properties
                mdl.elems(id).anm   = mdl.anm;
                mdl.elems(id).mat   = mdl.materials(p(1));
                mdl.elems(id).thk   = thickness(p(2));
                mdl.elems(id).shape = shape;
                
                % Always setup Gauss quadrature after setting shape because
                % element needs shape to define matrix for extrapolating
                % stresses from integration points to nodes
                mdl.elems(id).setGauss(gauss,intgrorder(p(3),1),intgrorder(p(3),2));
            end
        end
        
        %------------------------------------------------------------------
        function status = elementQuad4(this,fid,mdl,gauss,thickness,intgrorder)
            %                     4 +---------------+ 3
            %                       |               |
            %                       |               |
            %                       |     QUAD4     |
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
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nel,'number of elements Q4'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element Q4 ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Material ID, Thickness ID, Integration order ID, Connectivity
                [p,count] = fscanf(fid,'%d',7);
                if (~this.chkElemProp(mdl,id,p,thickness,intgrorder,count,7))
                    status = 0;
                    return;
                end
                
                % Create shape object
                shape = fem.Shape_Quad4([mdl.nodes(p(4));mdl.nodes(p(5));...
                                         mdl.nodes(p(6));mdl.nodes(p(7))]);
                
                % Store properties
                mdl.elems(id).anm   = mdl.anm;
                mdl.elems(id).mat   = mdl.materials(p(1));
                mdl.elems(id).thk   = thickness(p(2));
                mdl.elems(id).shape = shape;
                
                % Always setup Gauss quadrature after setting shape because
                % element needs shape to define matrix for extrapolating
                % stresses from integration points to nodes
                mdl.elems(id).setGauss(gauss,intgrorder(p(3),1),intgrorder(p(3),2));
            end
        end
        
        %------------------------------------------------------------------
        function status = elementTria6(this,fid,mdl,gauss,thickness,intgrorder)
            %                         3 +
            %                           |\
            %                           | \
            %                         6 +  + 5
            %                           |   \ TRIA6
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
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nel,'number of elements T6'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element T6 ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Material ID, Thickness ID, Integration order ID, Connectivity
                [p,count] = fscanf(fid,'%d',9);
                if (~this.chkElemProp(mdl,id,p,thickness,intgrorder,count,9))
                    status = 0;
                    return;
                end
                
                % Create shape object:
                % Node incidence is given in counter clockwise order around
                % element shape and is converted to corner first incidence
                shape = fem.Shape_Tria6([mdl.nodes(p(4));mdl.nodes(p(6));...
                                         mdl.nodes(p(8));mdl.nodes(p(5));...
                                         mdl.nodes(p(7));mdl.nodes(p(9))]);
                
                % Store properties
                mdl.elems(id).anm   = mdl.anm;
                mdl.elems(id).mat   = mdl.materials(p(1));
                mdl.elems(id).thk   = thickness(p(2));
                mdl.elems(id).shape = shape;
                
                % Always setup Gauss quadrature after setting shape because
                % element needs shape to define matrix for extrapolating
                % stresses from integration points to nodes
                mdl.elems(id).setGauss(gauss,intgrorder(p(3),1),intgrorder(p(3),2));
            end
        end
        
        %------------------------------------------------------------------
        function status = elementQuad8(this,fid,mdl,gauss,thickness,intgrorder)
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
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nel,'number of elements Q8'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element Q8 ID'))
                    status = 0;
                    return;
                end
                
                % Element properties:
                % Material ID, Thickness ID, Integration order ID, Connectivity
                [p,count] = fscanf(fid,'%d',11);
                if (~this.chkElemProp(mdl,id,p,thickness,intgrorder,count,11))
                    status = 0;
                    return;
                end
                
                % Create shape object:
                % Node incidence is given in counter clockwise order around
                % element shape and is converted to corner first incidence
                shape = fem.Shape_Quad8([mdl.nodes(p(4)); mdl.nodes(p(6));...
                                         mdl.nodes(p(8)); mdl.nodes(p(10));...
                                         mdl.nodes(p(5)); mdl.nodes(p(7));...
                                         mdl.nodes(p(9)); mdl.nodes(p(11))]);
                
                % Store properties
                mdl.elems(id).anm   = mdl.anm;
                mdl.elems(id).mat   = mdl.materials(p(1));
                mdl.elems(id).thk   = thickness(p(2));
                mdl.elems(id).shape = shape;
                
                % Always setup Gauss quadrature after setting shape because
                % element needs shape to define matrix for extrapolating
                % stresses from integration points to nodes
                mdl.elems(id).setGauss(gauss,intgrorder(p(3),1),intgrorder(p(3),2));
            end
        end
        
        %--------------------------------------------------------------------------
        function status = loadPoint(this,fid,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before point loads!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with point load
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nnp,'number of nodes with point load'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nnp,'node ID for point load specification'))
                    status = 0;
                    return;
                end
                
                % Point load values
                [load,count] = fscanf(fid,'%f',6);
                if (count ~= 6)
                    fprintf('Invalid point load of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).load = load;
            end
        end
        
        %------------------------------------------------------------------
        function status = nodePrescDispl(this,fid,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before prescribed displacements!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with prescribed displacements
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nnp,'number of nodes with prescribed displacements'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nnp,'node ID for prescribed displacement specification'))
                    status = 0;
                    return;
                end
                
                % Prescribed displacements values
                [disp,count] = fscanf(fid,'%f',6);
                if (count ~= 6)
                    fprintf('Invalid prescribed displacements of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).ebcDispl = disp;
            end
        end
        
        %--------------------------------------------------------------------------
        function status = nodePrescTemp(this,fid,mdl)
            status = 1;
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before prescribed temperatures!\n');
                status = 0;
                return;
            end
            
            % Total number of nodes with prescribed temperature
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,mdl.nnp,'number of nodes with prescribed temperature'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nnp,'node ID for prescribed temperature specification'))
                    status = 0;
                    return;
                end
                
                % Prescribed displacements values
                [temp,count] = fscanf(fid,'%f',1);
                if (count ~= 1)
                    fprintf('Invalid prescribed temperature of node %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).fixTemp = true;
                mdl.nodes(id).ebcTemp = temp;
            end
        end
        
        %--------------------------------------------------------------------------
        function status = loadLineUnif(this,fid,mdl)
            status = 1;
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before line loads!\n');
                status = 0;
                return;
            end
            
            % Total number of edges with line load
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of edges with line load'))
                status = 0;
                return;
            end
            
            a = zeros(n,1);
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element ID for line load specification'))
                    status = 0;
                    return;
                end
                
                % Line load information (node1,node2,loc_gbl,qx,qy,qz)
                [load,count] = fscanf(fid,'%f',6);
                if (count ~= 6)
                    fprintf('Invalid line load specification of element %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                % All loads are considered in global directions, so loc_gbl is ignored
                a(i) = id;
                j = sum(a(:)==id);
                mdl.elems(id).lineForce(j,:) = [load(1),load(2),load(4),load(5),load(6)];
            end
        end
        
        %--------------------------------------------------------------------------
        function status = loadDomainUnif(this,fid,mdl)
            status = 1;
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before domain loads!\n');
                status = 0;
                return;
            end
            
            % Total number of elements with doamain load
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of elements with domain load'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element ID for domain load specification'))
                    status = 0;
                    return;
                end
                
                % Domain load information (qx,qy,qz)
                [load,count] = fscanf(fid,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid domain load specification of element %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.elems(id).domainForce = load;
            end
        end
        
        %--------------------------------------------------------------------------
        function status = fluxLineUnif(this,fid,mdl)
            status = 1;
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before line fluxes!\n');
                status = 0;
                return;
            end
            
            % Total number of edges with line flux
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of edges with line flux'))
                status = 0;
                return;
            end
            
            a = zeros(n,1);
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element ID for line flux specification'))
                    status = 0;
                    return;
                end
                
                % Line flux information (node1,node2,q)
                [flux,count] = fscanf(fid,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid line flux specification of element %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                a(i) = id;
                j = sum(a(:)==id);
                mdl.elems(id).lineForce(j,:) = flux;
            end
        end
        
        %--------------------------------------------------------------------------
        function status = fluxDomainUnif(this,fid,mdl)
            status = 1;
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before domain heat fluxes!\n');
                status = 0;
                return;
            end
            
            % Total number of elements with domain heat flux
            n = fscanf(fid,'%d',1);
            if (~this.chkInt(n,inf,'number of elements with domain heat flux'))
                status = 0;
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fid,'%d',1);
                if (~this.chkInt(id,mdl.nel,'element ID for domain flux specification'))
                    status = 0;
                    return;
                end
                
                % Domain flux value
                [flux,count] = fscanf(fid,'%f',1);
                if (count ~= 1)
                    fprintf('Invalid domain flux specification of element %d\n',id);
                    status = 0;
                    return;
                end
                
                % Store data
                mdl.elems(id).domainForce = flux;
            end
        end
        
        %% Functions for checking input data
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
        function status = chkElemProp(this,mdl,id,prop,thks,intord,count,num)
            status = 0;
            nodesIDs = prop(4:end);
            if (count ~= num)
                fprintf('Invalid properties of element %d\n',id);
            elseif (~this.chkInt(prop(1),mdl.nmat,'material ID for element definition'))
            elseif (~this.chkInt(prop(2),size(thks,1),'thickness ID for element definition'))
            elseif (~this.chkInt(prop(3),size(intord,1),'integration order ID for element definition'))
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