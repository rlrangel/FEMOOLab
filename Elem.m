%% Elem class (Element class)
classdef Elem < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 0;                  % identification number
        elem_type = 0;           % flag for element type (T3=0, Q4=1, T6=2, Q8=3)
        edge_type = 0;           % type of finite element edge (linear=0, quadratic=1)
        gauss_type = 0;          % type of Gauss quadrature (line=0, trinagular=1, quadrilateral=2)
        gauss_order = 0;         % order of Gauss quadrature for stiffness matrix computation
        nen = 0;                 % number of nodes
        nedgen = 0;              % number of edge nodes
        anm = [];                % handle to an object of the Anm class
        nodes = [];              % vector of handles to objects of the Node class
        material = [];           % handle to an object of the Material class
        thk = 0;                 % element thickness
        gle = [];                % gather vector (stores element d.o.f. eqn. numbers)
        edgeLoad = [];           % 5 column matrix with edge uniform distrib. load components [node_i node_j qx qy qz]
        areaLoad = [];           % 3 position vector with area uniform distrib. load components [px, py, pz]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem(type,edge,gauss,nen,nedgen)
            elem.elem_type = type;
            elem.edge_type = edge;
            elem.gauss_type = gauss;
            elem.nen = nen;
            elem.nedgen = nedgen;
        end
    end
    
    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Evaluates shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  N:    shape function matrix evaluated at given position
        N = shpNmatrix(elem,r,s)
        
        %------------------------------------------------------------------
        % Evaluates derivatives of shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  GradNpar: Shape function gradient matrix evaluated at given position
        GradNpar = shpGradNmatrix(elem,r,s)
        
        %------------------------------------------------------------------
        % Get Gauss points coordinates in the parametric space of element and
        % corresponding weights for a given quadrature type and order.
        % Output arguments:
        %  ngp: number of Gauss points
        %  w: vector of Gauss point weights
        %  gp: vector of Gauss point parametric coordinates
        [ngp,w,gp] = gauss(elem)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes element stiffness matrix.
        % Output:
        %  ke: element stiffness matrix
        function ke = stiffMtx(elem)
            % Initialize element stiffness matrix
            ke = zeros(elem.nen * elem.anm.ndof, elem.nen * elem.anm.ndof);
            
            % Generate material constituive matrix
            E = elem.anm.EMtx(elem);
            
            % Generate cartesian nodal coordinates matrix
            X = elem.anm.elemNodalCoord(elem);
            
            % Get gauss points and weights
            [ngp,w,gp] = elem.gauss();
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Get gauss point parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Calculate matrix of shape functions derivatives in parametric coordinates
                GradNpar = elem.shpGradNmatrix(r,s);
                
                % Compute Jacobian matrix evaluated at current Gauss point
                J = GradNpar * X;
                
                % Calculate matrix of shape functions derivatives in cartesian coordinates
                GradNcar = J \ GradNpar;
                
                % Compute B matrix evaluated at current Gauss point
                B = elem.anm.BMtx(elem,GradNcar,X,r,s);
                
                % Accumulate gauss point contribution to element stiffness matrix
                ke = ke + w(i) * elem.thk * det(J) * B' * E * B;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element equivalent nodal load vector for area loads.
        % Output:
        %  fe: equivalent nodal load vector
        function fe = areaENL(elem)
            % Initialize element equiv. nodal load vector
            fe = zeros(elem.nen * elem.anm.ndof, 1);
            
            % Generate cartesian nodal coordinates matrix
            X = elem.anm.elemNodalCoord(elem);
            
            % Get gauss points and weights
            [ngp,w,gp] = elem.gauss();
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Get gauss point parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Compute shape functions matrix evaluated at this gauss point
                N = elem.shpNmatrix(r,s);
                
                % Calculate matrix of shape functions derivatives in parametric coordinates
                GradNpar = elem.shpGradNmatrix(r,s);
                
                % Compute Jacobian matrix evaluated at current Gauss point
                J = GradNpar * X;
                detJ = det(J);
                
                % Accumulate gauss point contribution
                il = 0;
                for j = 1:elem.nen                     
                    for k = 1:elem.anm.ndof
                        il = il + 1;
                        fe(il) = fe(il) + w(i) * detJ * N(j) * elem.areaLoad(k);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Computes element equivalent nodal load vector for edge loads.
        % Input arguments:
        %  l: line of edge load matrix of current element
        % Output:
        %  fe: equivalent nodal load vector
        function fe = edgeENL(elem,l)
            % Initialize element equiv. nodal load vector
            fe = zeros(elem.nedgen * elem.anm.ndof, 1);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialize element equiv. nodal load vector
            fe = zeros(elem.nen * elem.anm.ndof, 1);
            
            % Generate cartesian nodal coordinates matrix
            X = elem.anm.elemNodalCoord(elem);
            
            % Get gauss points and weights
            [ngp,w,gp] = elem.gauss();
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Get gauss point parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Compute shape functions matrix evaluated at this gauss point
                N = elem.shpNmatrix(r,s);
                
                % Calculate matrix of shape functions derivatives in parametric coordinates
                GradNpar = elem.shpGradNmatrix(r,s);
                
                % Compute Jacobian matrix evaluated at current Gauss point
                J = GradNpar * X;
                detJ = det(J);
                
                % Accumulate gauss point contribution
                il = 0;
                for j = 1:elem.nen                     
                    for k = 1:elem.anm.ndof
                        il = il + 1;
                        fe(il) = fe(il) + w(i) * detJ * N(j) * elem.areaLoad(k);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of an Elem object.
        function clean(elem)
            elem.id = 0;
            elem.elem_type = 0;
            elem.edge_type = 0;
            elem.gauss_type = 0;
            elem.gauss_order = 0;
            elem.nen = 0;
            elem.nedgen = 0;
            elem.anm = [];
            elem.nodes = [];
            elem.material = [];
            elem.thk = 0;
            elem.gle = [];
            elem.edgeLoad = [];
            elem.areaLoad = [];
        end
    end
end