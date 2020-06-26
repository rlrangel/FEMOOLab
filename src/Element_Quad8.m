%% Element_Quad8 Class (Serendipity Quadratic Quadrilateral Element)
%
% This is a sub-class in the StAnOOP program that implements abstract 
% methods declared in super-class Element to deal with 8-noded quadrilateral
% (serendipity quadratic quadrilateral) elements:
%                                   7
%                         4 +-------+-------+ 3
%                           |       s       |
%                           |       ^       |
%                           |       |       |
%                         8 +       ----> r + 6
%                           |               |
%                           |     QUAD8     |
%                           |               |
%                         1 +-------+-------+ 2
%                                   5
classdef Element_Quad8 < Element
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Element_Quad8(id,anm,thk,mat,gstiff_order,gstress_order,nodes)
            c = Constants();
            elem = elem@Element(c.QUAD8,8);
            
            if (nargin > 0)
                elem.id            = id;
                elem.anm           = anm;
                elem.thk           = thk;
                elem.mat           = mat;
                elem.nodes         = nodes;
                elem.edge_type     = c.LINE3;
                elem.nedgen        = 3;
                elem.gauss_type    = c.QUADRATURE_QUAD;
                elem.gstiff_order  = gstiff_order;
                elem.gstress_order = gstress_order;
                
                % Cartesian nodal coordiantes matrix [X Y]
                elem.carCoord =...
                [ nodes(1).coord(1)   nodes(1).coord(2);
                  nodes(2).coord(1)   nodes(2).coord(2);
                  nodes(3).coord(1)   nodes(3).coord(2);
                  nodes(4).coord(1)   nodes(4).coord(2);
                  nodes(5).coord(1)   nodes(5).coord(2);
                  nodes(6).coord(1)   nodes(6).coord(2);
                  nodes(7).coord(1)   nodes(7).coord(2);
                  nodes(8).coord(1)   nodes(8).coord(2) ];
                
                % Parametric nodal coordinates matrix [r s]
                elem.parCoord = [ -1 -1;
                                   1 -1;
                                   1  1;
                                  -1  1;
                                   0 -1;
                                   1  0;
                                   0  1;
                                  -1  0 ];
                
                % Number of Gauss points for stress evaluation and
                % Transformation matrix of Gauss-to-node results
                [ngp,TGN] = elem.TGNmtx();
                elem.gstress_npts = ngp;
                elem.TGN = TGN;
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Element
    methods
        %------------------------------------------------------------------
        % Get parametric coordinates and weights of Gauss points for a
        % given quadrature order.
        % Output:
        %  ngp: number of Gauss points
        %  w:   vector of Gauss points weights
        %  gp:  array of Gauss points parametric coordinates
        function [ngp,w,gp] = gauss(~,order)  
            % Number of Gauss points
            ngp = order * order;
            
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
            npnts = 0;
            for i = 1:order
                for j = 1:order
                    npnts = npnts + 1;
                    w(npnts)    = qwgt(j,order) * qwgt(i,order);
                    gp(1,npnts) = qpts(i,order);
                    gp(2,npnts) = qpts(j,order);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of shape functions at a given position in
        % parametric coordinates.
        function N = Nmtx(elem,r,s)
            N = zeros(1,elem.nen);
            
            N(5) = 0.50*(1-r*r)*(1-s);
            N(6) = 0.50*(1+r)*(1-s*s);
            N(7) = 0.50*(1-r*r)*(1+s);
            N(8) = 0.50*(1-r)*(1-s*s);
            N(1) = 0.25*(1-r)*(1-s) - 0.50*N(8) - 0.50*N(5);
            N(2) = 0.25*(1+r)*(1-s) - 0.50*N(5) - 0.50*N(6);
            N(3) = 0.25*(1+r)*(1+s) - 0.50*N(6) - 0.50*N(7);
            N(4) = 0.25*(1-r)*(1+s) - 0.50*N(7) - 0.50*N(8);
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge shape functions at a given position in
        % parametric coordinate.
        function N = NmtxEdge(elem,r)
            N = zeros(1,elem.nedgen);
            
            N(1) = 0.5*(1-r) - 0.5*N(3);
            N(2) = 0.5*(1+r) - 0.5*N(3);
            N(3) = 1-r*r;
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of shape functions derivatives w.r.t. parametric
        % coordinates at a given position.
        function GradNpar = gradNmtx(elem,r,s)
            GradNpar = zeros(2,elem.nen);
            
            GradNpar(1,1) = (2*r - 2*r*s - s*s   + s) / 4;
            GradNpar(2,1) = (2*s - r*r   - 2*r*s + r) / 4;
            GradNpar(1,2) = (2*r - 2*r*s + s*s   - s) / 4;
            GradNpar(2,2) = (2*s - r*r   + 2*r*s - r) / 4;
            GradNpar(1,3) = (2*r + 2*r*s + s*s   + s) / 4;
            GradNpar(2,3) = (2*s + r*r   + 2*r*s + r) / 4;
            GradNpar(1,4) = (2*r + 2*r*s - s*s   - s) / 4;
            GradNpar(2,4) = (2*s + r*r   - 2*r*s - r) / 4;
            GradNpar(1,5) = r*s - r;
            GradNpar(2,5) = (r*r - 1) / 2;
            GradNpar(1,6) = (1 - s*s) / 2;
            GradNpar(2,6) = -s - r*s;
            GradNpar(1,7) = -r - r*s;
            GradNpar(2,7) = (1 - r*r) / 2;
            GradNpar(1,8) = (s*s - 1) / 2;
            GradNpar(2,8) = r*s - s;
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of edge shape functions derivatives w.r.t. parametric
        % coordinates at a given position.
        function GradNpar = gradNmtxEdge(elem,r)
            GradNpar = zeros(1,elem.nedgen);
            
            GradNpar(1,1) = -0.5 + r;
            GradNpar(1,2) =  0.5 + r;
            GradNpar(1,3) = -2*r;
        end
        
        %------------------------------------------------------------------
        % Compute equivalent nodal load (ENL) vector for line loads
        % distributed over the edge of elements.
        function f = lineEquivLoadVct(elem)
            ndof   = elem.anm.ndof;
            nen    = elem.nen;
            nedgen = elem.nedgen;
            
            % Initialize element ENL vector
            f = zeros(nen*ndof,1);
            
            % Loop over element line loads
            for q = 1:size(elem.lineLoad,1)
                % Initialize edge ENL vector
                fline = zeros(nedgen*ndof,1);
                
                % Global IDs of edge initial and final nodes
                n1 = elem.lineLoad(q,1);
                n3 = elem.lineLoad(q,2);
                
                % Local IDs of edge nodes
                for i = 1:nen
                    % Vertice node IDs
                    node = elem.nodes(i).id;
                    if node == n1
                        n1 = i;
                    elseif node == n3
                        n3 = i;
                    end
                    
                    % Middle edge node ID
                    if (n1 == 1 || n3 == 1) && (n1 == 2 || n3 == 2)
                        n2 = 5;
                    elseif (n1 == 2 || n3 == 2) && (n1 == 3 || n3 == 3)
                        n2 = 6;
                    elseif (n1 == 3 || n3 == 3) && (n1 == 4 || n3 == 4)
                        n2 = 7;
                    elseif (n1 == 1 || n3 == 1) && (n1 == 4 || n3 == 4)
                        n2 = 8;
                    end
                end
                
                % Edge nodes coordinates
                X = elem.carCoord([n1,n2,n3],:);
                
                % Load components
                p = elem.lineLoad(q,3:5);
                
                % Gauss points and weights
                [ngp,w,gp] = elem.gaussEdge(elem.gstiff_order);
                
                % Loop over edge Gauss integration points
                for i = 1:ngp
                    % Parametric coordinates
                    r = gp(1,i);
                    
                    % Edge shape functions matrix
                    N = elem.NmtxEdge(r);
                    
                    % Matrix of edge shape functions derivatives w.r.t. parametric coordinates
                    GradNpar = elem.gradNmtxEdge(r);
                    
                    % Jacobian matrix
                    J = GradNpar * X;
                    detJ = sqrt(J(1)*J(1) + J(2)*J(2));
                    
                    % Accumulate gauss point contributions
                    m = 0;
                    for j = 1:nedgen
                        for k = 1:ndof
                            m = m + 1;
                            fline(m) = fline(m) + w(i) * detJ * N(j) * p(k);
                        end
                    end
                end
                
                % Edge gather vector (stores local d.o.f.'s numbers)
                gledge = [ndof*n1-1 ndof*n1 ndof*n2-1 ndof*n2 ndof*n3-1 ndof*n3];
                
                % Assemble edge ENL vetor to element ENL vector
                f(gledge) = f(gledge) + fline;
            end
        end
    end
end