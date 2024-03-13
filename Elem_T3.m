%% Elem_T3 class (Linear Triangle Element)
%                           s
%                           ^
%                           |
%
%                         3 +
%                           |\
%                           | \
%                           |  \
%                           |   \ TRIA3
%                           |    \
%                         1 +-----+ 2  ----> r
%
classdef Elem_T3 < Elem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem_T3()
            include_gblrefs;
            elem = elem@Elem(TRIA3,LINE2,TRIA_QUADRATURE,3,2);
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
            
            N(1) = 1 - r - s;
            N(2) = r;
            N(3) = s;
        end
        
        %------------------------------------------------------------------
        % Evaluates derivatives of shape functions in parametric coordinates.
        % Input arguments:
        %  r,s:  parametric coordinate values
        % Output arguments:
        %  GradNpar: Shape function gradient matrix evaluated at given position
        function GradNpar = shpGradNmatrix(~,~,~)
            
            GradNpar = [ -1  1  0;
                         -1  0  1 ];
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