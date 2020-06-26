%% Element_Tri6 Class (Quadratic Triangular Element)
%
% This is a sub-class in the StAnOOP program that implements abstract 
% methods declared in super-class Element to deal with 6-noded trianlge
% (quadratic triangular) elements:
%                           s
%                           ^
%                           |
%
%                         3 +
%                           |\
%                           | \
%                         6 +  + 5
%                           |   \ TRI6
%                           |    \
%                         1 +--+--+ 2  ----> r
%                              4
classdef Element_Tri6 < Element
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Element_Tri6(id,anm,thk,mat,gstiff_order,gstress_order,nodes)
            c = Constants();
            elem = elem@Element(c.TRI6,6);
            
            if (nargin > 0)
                elem.id            = id;
                elem.anm           = anm;
                elem.thk           = thk;
                elem.mat           = mat;
                elem.nodes         = nodes;
                elem.edge_type     = c.LINE3;
                elem.nedgen        = 3;
                elem.gauss_type    = c.QUADRATURE_TRIA;
                elem.gstiff_order  = gstiff_order;
                elem.gstress_order = gstress_order;
                
                % Cartesian nodal coordiantes matrix [X Y]
                elem.carCoord =...
                [ nodes(1).coord(1)   nodes(1).coord(2);
                  nodes(2).coord(1)   nodes(2).coord(2);
                  nodes(3).coord(1)   nodes(3).coord(2);
                  nodes(4).coord(1)   nodes(4).coord(2);
                  nodes(5).coord(1)   nodes(5).coord(2);
                  nodes(6).coord(1)   nodes(6).coord(2) ];
                
                % Parametric nodal coordinates matrix [r s]
                elem.parCoord = [ 0.0  0.0;
                                  1.0  0.0;
                                  0.0  1.0;
                                  0.5  0.0;
                                  0.5  0.5;
                                  0.0  0.5 ];
                
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
            if (order == 1)
                ngp = 1;
                
                w = 0.5;
                
                gp = [ 0.333333333333
                       0.333333333333 ];
                
            elseif (order == 2 || order == 3)
                ngp = 3;
                
                w  = [ 0.166666666666 0.166666666666 0.166666666666 ];
                
                gp = [ 0.166666666666 0.666666666666 0.166666666666
                       0.166666666666 0.166666666666 0.666666666666 ];
                
            elseif (order == 4)
                ngp = 4;
                
                w = [ -0.281250000000 0.260416666666 0.260416666666 0.260416666666 ];
                
                gp = [ 0.333333333333 0.200000000000 0.600000000000 0.200000000000
                       0.333333333333 0.200000000000 0.200000000000 0.600000000000 ];
            end
        end
        
        %------------------------------------------------------------------
        % Evaluate matrix of shape functions at a given position in
        % parametric coordinates.
        function N = Nmtx(elem,r,s)
            N = zeros(1,elem.nen);
            
            N(1) = 1 - 3*r - 3*s + 4*r*s + 2*r*r + 2*s*s;
            N(2) = 2*r*r - r;
            N(3) = 2*s*s - s;
            N(4) = 4*r - 4*r*r - 4*r*s;
            N(5) = 4*r*s;
            N(6) = 4*s - 4*r*s - 4*s*s;
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
            
            GradNpar(1,1) =  4*r + 4*s - 3;
            GradNpar(2,1) =  4*r + 4*s - 3;
            GradNpar(1,2) =  4*r - 1;
            GradNpar(2,2) =  0;
            GradNpar(1,3) =  0;
            GradNpar(2,3) =  4*s - 1;
            GradNpar(1,4) =  4 - 8*r - 4*s;
            GradNpar(2,4) = -4*r;
            GradNpar(1,5) =  4*s;
            GradNpar(2,5) =  4*r;
            GradNpar(1,6) = -4*s;
            GradNpar(2,6) =  4 - 4*r - 8*s;
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