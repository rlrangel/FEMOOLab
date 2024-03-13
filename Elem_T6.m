%% Elem_T6 class (Quadratic Triangle Element)
%                           s
%                           ^
%                           |
%
%                         3 +
%                           |\
%                           | \
%                         6 +  + 5
%                           |   \ TRIA6
%                           |    \
%                         1 +--+--+ 2  ----> r
%                              4
classdef Elem_T6 < Elem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem_T6()
            include_gblrefs;
            elem = elem@Elem(TRIA6,LINE3,TRIA_QUADRATURE,6,3);
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
            
            N(1) = 1 - 3*r - 3*s + 4*r*s + 2*r*r + 2*s*s;
            N(2) = 2*r*r - r;
            N(3) = 2*s*s - s;
            N(4) = 4*r - 4*r*r - 4*r*s;
            N(5) = 4*r*s;
            N(6) = 4*s - 4*r*s - 4*s*s;
        end
        
        %------------------------------------------------------------------
        % Evaluates derivatives of shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  GradNpar: Shape function gradient matrix evaluated at given position
        function GradNpar = shpGradNmatrix(elem,r,s)
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
        % Get Gauss points coordinates in the parametric space of element and
        % corresponding weights for a given quadrature type and order.
        % Output arguments:
        %  ngp: number of Gauss points
        %  w: vector of Gauss point weights
        %  gp: vector of Gauss point parametric coordinates
        function [ngp,w,gp] = gauss(elem)
            % Assemble arrays of weights and gauss point coordinates
            if elem.gauss_order == 1
                % Number of Gauss points
                ngp = 1;
                
                % Assemble array of weights
                w = 0.500000000000;
                
                % Assemble array of gauss point coordinates
                gp = [ 0.333333333333
                       0.333333333333 ];
                
            elseif (elem.gauss_order == 2) || (elem.gauss_order == 3)
                % Number of Gauss points
                ngp = 3;
                
                % Assemble array of weights
                w = [ 0.166666666666  0.166666666666 0.166666666666 ];
                
                % Assemble array of gauss point coordinates
                gp = [ 0.166666666666 0.666666666666 0.166666666666
                       0.166666666666 0.166666666666 0.666666666666 ];
                
            elseif elem.gauss_order == 4
                % Number of Gauss points
                ngp = 4;
                
                % Assemble array of weights
                w = [ -0.281250000000  0.260416666666 0.260416666666 0.260416666666 ];
                
                % Assemble array of gauss point coordinates
                gp = [ 0.333333333333 0.200000000000 0.600000000000 0.200000000000
                       0.333333333333 0.200000000000 0.200000000000 0.600000000000 ];
                
            end
        end
    end
end