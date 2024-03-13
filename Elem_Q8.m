%% Elem_Q8 class (Serendipity Quadratic Quadrilateral Element)
%
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
classdef Elem_Q8 < Elem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem_Q8()
            include_gblrefs;
            elem = elem@Elem(QUAD8,LINE3,QUAD_QUADRATURE,8,3);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <elem.html *Elem*>.
    methods
        %------------------------------------------------------------------
        % Evaluates shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  N:    shape function matrix evaluated at given position
        function N = shpNmatrix(elem,r,s)
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
        % Evaluates derivatives of shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  GradNpar: Shape function gradient matrix evaluated at given position
        function GradNpar = shpGradNmatrix(elem,r,s)
            GradNpar = zeros(2,elem.nen);
            r2 = r*r;
            s2 = s*s;
            rs = r*s;
            
            GradNpar(1,1) = (2*r - 2*rs - s2 + s) / 4;
            GradNpar(2,1) = (2*s - r2 - 2*rs + r) / 4;
            GradNpar(1,2) = (2*r - 2*rs + s2 - s) / 4;
            GradNpar(2,2) = (2*s - r2 + 2*rs - r) / 4;
            GradNpar(1,3) = (2*r + 2*rs + s2 + s) / 4;
            GradNpar(2,3) = (2*s + r2 + 2*rs + r) / 4;
            GradNpar(1,4) = (2*r + 2*rs - s2 - s) / 4;
            GradNpar(2,4) = (2*s + r2 - 2*rs - r) / 4;
            GradNpar(1,5) = rs - r;
            GradNpar(2,5) = (r2 - 1) / 2;
            GradNpar(1,6) = (1 - s2) / 2;
            GradNpar(2,6) = -s - rs;
            GradNpar(1,7) = -r - rs;
            GradNpar(2,7) = (1 - r2) / 2;
            GradNpar(1,8) = (s2 - 1) / 2;
            GradNpar(2,8) = rs - s;
        end
        
        %------------------------------------------------------------------
        %% Get Gauss points coordinates in the parametric space of element and
        % corresponding weights for a given quadrature type and order.
        % Output arguments:
        %  ngp: number of Gauss points
        %  w: vector of Gauss point weights
        %  gp: vector of Gauss point parametric coordinates
        function [ngp,w,gp] = gauss(elem)
            % Number of Gauss points
            order = elem.gauss_order;
            ngp = order * order;
            
            % Initialize arrays of Gauss point parametric coordenates and weights
            gp = zeros(2,ngp);
            w = zeros(1,ngp);
            
            % Matrix of weights
            qwgt = [ 2.00000000000000   1.000000000000000  0.555555555555555  0.347854845137454
                     0.00000000000000   1.000000000000000  0.888888888888888  0.652145154862546
                     0.00000000000000   0.000000000000000  0.555555555555555  0.652145154862546
                     0.00000000000000   0.000000000000000  0.000000000000000  0.347854845137454 ];
            
            % Matrix of Gauss points coordinates
            qpts = [ 0.000000000000000 -0.577350269189626 -0.774596669241483 -0.861136311594053
                     0.000000000000000  0.577350269189626  0.000000000000000 -0.339981043584856
                     0.000000000000000  0.000000000000000  0.774596669241483  0.339981043584856
                     0.000000000000000  0.000000000000000  0.000000000000000  0.861136311594053 ];
            
            % Assemble arrays of gauss point coordinates and weights
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
    end
end