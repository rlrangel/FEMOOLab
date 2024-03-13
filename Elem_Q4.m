%% Elem_Q4 class (Bilinear Quadrilateral Element)
%
%
%                         4 +---------------+ 3
%                           |       s       |
%                           |       ^       |
%                           |       |       |
%                           |       ----> r |
%                           |               |
%                           |     QUAD4     |
%                           |               |
%                         1 +---------------+ 2
%
classdef Elem_Q4 < Elem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem_Q4()
            include_gblrefs;
            elem = elem@Elem(QUAD4,LINE2,QUAD_QUADRATURE,4,2);
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
            
            N(1) = 0.25*(1-r)*(1-s);
            N(2) = 0.25*(1+r)*(1-s);
            N(3) = 0.25*(1+r)*(1+s);
            N(4) = 0.25*(1-r)*(1+s);
        end
        
        %------------------------------------------------------------------
        % Evaluates derivatives of shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  GradNpar: Shape function gradient matrix evaluated at given position
        function GradNpar = shpGradNmatrix(elem,r,s)
            GradNpar = zeros(2,elem.nen);
            
            GradNpar(1,1) = -0.25*(1 - s);
            GradNpar(2,1) = -0.25*(1 - r);
            GradNpar(1,2) =  0.25*(1 - s);
            GradNpar(2,2) = -0.25*(1 + r);
            GradNpar(1,3) =  0.25*(1 + s);
            GradNpar(2,3) =  0.25*(1 + r);
            GradNpar(1,4) = -0.25*(1 + s);
            GradNpar(2,4) =  0.25*(1 - r);
        end
        
        %------------------------------------------------------------------
        % Get Gauss points coordinates in the parametric space of element and
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