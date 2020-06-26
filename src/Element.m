%% Element Class
%
% This is an abstract super-class that generically specifies an element 
% in the StAnOOP program.
%
% Essentially, this super-class declares abstract methods that define the
% particular behavior of an element. These abstract methods are the
% functions that should be implemented in a derived sub-class that deals
% with specific types of elements.
%
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General:
        id   = 0;            % identification number
        type = 0;            % flag for element type 
        anm  = [];           % object of the Anm (Analysis Model) class
        gle  = [];           % gather vector (stores element global d.o.f. numbers)
        
        % Physical attributes
        thk = 0;             % thickness
        mat = [];            % object of the Material class
        
        % Geometry
        nen      = 0;        % number of nodes
        nodes    = [];       % vector of objects of the Node class
        carCoord = [];       % matrix of cartesian nodal coordinates
        parCoord = [];       % matrix of parametric nodal coordinates
        
        % Edges
        edge_type = 0;       % flag for edge type (linear=0, quadratic=1)
        nedgen    = 0;       % number of edge nodes
        
        % Gauss points
        gauss_type    = 0;   % flag for type of Gauss quadrature (line, triang., or quad.)
        gstiff_order  = 0;   % order of Gauss quadrature for stiffness matrix computation
        gstress_order = 0;   % order of Gauss quadrature for stress computation
        gstress_npts  = 0;   % number of gauss points for stress computation
        TGN           = [];  % transformation matrix of gauss results to node results
        
        % Loads
        lineLoad = [];       % matrix of uniform line loads [node1,node1,qx,qy,qz]
        areaLoad = [];       % vector of uniform area loads [px,py,pz]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Element(type,nen)
            elem.type = type;
            elem.nen  = nen;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Get parametric coordinates and weights of Gauss points for a
        % given quadrature order.
        [ngp,w,gp] = gauss(elem,order);
        
        %------------------------------------------------------------------
        % Evaluate matrix of shape functions at a given position in
        % parametric coordinates.
        N = Nmtx(elem,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge shape functions at a given position in
        % parametric coordinate.
        N = NmtxEdge(elem,r);
        
        %------------------------------------------------------------------
        % Evaluate matrix of shape functions derivatives w.r.t. parametric
        % coordinates at a given position.
        GradNpar = gradNmtx(elem,r,s);
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge shape functions derivatives w.r.t. parametric
        % coordinates at a given position.
        GradNpar = gradNmtxEdge(elem,r);
        
        %------------------------------------------------------------------
        % Compute equivalent nodal load (ENL) vector for line loads
        % distributed over the edge of elements.
        f = lineEquivLoadVct(elem);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute transformation matrix of gauss-to-node results .
        % Refs.:
        % -Hinton & Campbell, "Local and Global Smoothing of Discontinous
        % Finite Element Functions using a Least Squares Method",
        % Int. J. Num. Meth. Engng., Vol. 8, pp. 461-480, 1974.
        % -Burnett, D.S., "Finite Element Analysis - From Concepts to Applications",
        % Addison-Wesley, 1987.
        % -Martha, L.F., "Notas de Aula do Curso CIV 2118 - Metodo dos Elementos
        % Finitos", 1994.
        function [ngp,TGN] = TGNmtx(elem)
            % Get parametric coordinates of Gauss points
            [ngp,~,gp] = elem.gauss(elem.gstress_order);
            
            if (elem.gstress_order == 1)
                TGN = ones(elem.nen,1);
            else
                % Compute S matrix that defines the coefficients of the
                % smoothing stress plane that fits the Gauss point stress
                % values in a least square sence
                P = zeros(3,3);
                P(1,1) = ngp;
                P(1,2) = sum(gp(1,:));
                P(1,3) = sum(gp(2,:));
                P(2,1) = P(1,2);
                P(2,2) = gp(1,:) * gp(1,:)';
                P(2,3) = gp(1,:) * gp(2,:)';
                P(3,1) = P(1,3);
                P(3,2) = P(2,3);
                P(3,3) = gp(2,:) * gp(2,:)';
                
                Q = ones(3,ngp);
                Q(2,:) = gp(1,:);
                Q(3,:) = gp(2,:);
                
                S = P\Q;
                
                % Compute nodal stress evaluation matrix, which is obtained
                % using nodal parametric coordinates in the smoothing stress
                % plane equation
                E = ones(elem.nen,3);
                E(:,2) = elem.parCoord(:,1);
                E(:,3) = elem.parCoord(:,2);
                
                % Compute transformation matrix
                TGN = E * S;
            end
        end
        
        %------------------------------------------------------------------
        % Compute elastic stiffness matrix.
        function k = elastStiffMtx(elem)
            % Initialize element matrix
            k = zeros(elem.nen*elem.anm.ndof, elem.nen*elem.anm.ndof);
            
            % Cartesian coordinates matrix
            X = elem.carCoord;
            
            % Material constitutive matrix
            C = elem.anm.Cmtx(elem);
            
            % Gauss points and weights
            [ngp,w,gp] = elem.gauss(elem.gstiff_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Matrix of shape functions derivatives w.r.t. parametric coordinates
                GradNpar = elem.gradNmtx(r,s);
                
                % Jacobian matrix
                J = GradNpar * X;
                
                % Matrix of shape functions derivatives w.r.t. cartesian coordinates
                GradNcar = J \ GradNpar;
                
                % Strain-displacemen matrix 
                B = elem.anm.Bmtx(elem,GradNcar,r,s);
                
                % Accumulate Gauss point contributions
                k = k + w(i) * elem.thk * det(J) * B' * C * B;
            end
        end
        
        %------------------------------------------------------------------
        % Compute equivalent nodal load (ENL) vector for area loads.
        function f = areaEquivLoadVct(elem)
            % Initialize element ENL vector
            f = zeros(elem.nen*elem.anm.ndof,1);
            
            % Cartesian coordinates matrix
            X = elem.carCoord;
            
            % Gauss points and weights
            [ngp,w,gp] = elem.gauss(elem.gstiff_order);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = elem.Nmtx(r,s);
                
                % Matrix of shape functions derivatives w.r.t. parametric coordinates
                GradNpar = elem.gradNmtx(r,s);
                
                % Jacobian matrix
                J = GradNpar * X;
                detJ = det(J);
                
                % Accumulate gauss point contributions
                m = 0;
                for j = 1:elem.nen
                    for k = 1:elem.anm.ndof
                        m = m + 1;
                        f(m) = f(m) + w(i) * detJ * N(j) * elem.areaLoad(k);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Get parametric coordinates and weights of edge Gauss points for a
        % linear quadrature and given order.
        % Input:
        %  order: quadrature order
        % Output:
        %  ngp: number of Gauss points
        %  w:   vector of Gauss points weights
        %  gp:  array of Gauss points parametric coordinates
        function [ngp,w,gp] = gaussEdge(~,order)
            % Number of Gauss points
            ngp = order;
            
            % Initialize arrays of weights and parametric coordinates
            w  = zeros(1,ngp);
            gp = zeros(2,ngp);
            
            % Matrix of weights
            qwgt = [ 2.00000000000000   1.000000000000000  0.555555555555555  0.347854845137454
                     0.00000000000000   1.000000000000000  0.888888888888888  0.652145154862546
                     0.00000000000000   0.000000000000000  0.555555555555555  0.652145154862546
                     0.00000000000000   0.000000000000000  0.000000000000000  0.347854845137454 ];
            
            % Matrix of parametric coordinates
            qpts = [ 0.000000000000000 -0.577350269189626 -0.774596669241483 -0.861136311594053
                     0.000000000000000  0.577350269189626  0.000000000000000 -0.339981043584856
                     0.000000000000000  0.000000000000000  0.774596669241483  0.339981043584856
                     0.000000000000000  0.000000000000000  0.000000000000000  0.861136311594053 ];
            
            % Assemble arrays of weights and parametric coordinates
            for i = 1:order
                w(i)    = qwgt(i,order);
                gp(1,i) = qpts(i,order);
                gp(2,i) = 0.0;
            end
        end
        
        %------------------------------------------------------------------
        % Compute stress components at gauss points and gauss point cartesian
        % coordinates for a given element.
        % Input:
        %  d: solution at element nodes
        % Output:
        %  ngp:  number of Gauss points for stress evaluation
        %  str:  stress components (sx,sy,txy) at each Gauss point
        %  gpc:  Gauss point cartesian coordinates array
        function [ngp,str,gpc] = gaussStress(elem,d)
            % Get Gauss points and weights
            [ngp,~,gp] = elem.gauss(elem.gstress_order);
            
            % Initialize stress component matrix and Gauss point coordinates
            str = zeros(3,ngp);
            gpc = zeros(2,ngp);
            
            % Cartesian coordinates
            X = elem.carCoord;
            
            % Material constituive matrix
            C = elem.anm.Cmtx(elem);
            
            % Loop over Gauss integration points
            for i = 1:ngp
                % Parametric coordinates
                r = gp(1,i);
                s = gp(2,i);
                
                % Shape functions matrix
                N = elem.Nmtx(r,s);
                
                % Matrix of shape functions derivatives w.r.t. parametric coordinates
                GradNpar = elem.gradNmtx(r,s);
                
                % Jacobian matrix
                J = GradNpar * X;
                
                % Matrix of shape functions derivatives w.r.t. cartesian coordinates
                GradNcar = J \ GradNpar;
                
                % Strain-displacemen matrix 
                B = elem.anm.Bmtx(elem,GradNcar,r,s);
                
                % Gauss point stress components and cartesian coordinates
                str(:,i) = elem.anm.pointStress(C,B,d);
                gpc(:,i) = N * X;
            end
        end
    end
end