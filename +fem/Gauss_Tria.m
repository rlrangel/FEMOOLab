%% Gauss_Tria Class (Gauss Triangular Quadrature)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <gauss.html Gauss: integration quadrature super-class>
% to deal with triangular element shapes.
%
% Eg.: Quadrature of order 2 or 3:
%                              s
%                              ^
%                              |
%                              +
%                              |\
%                              | \
%                              |  \
%                              | * \
%                              | 3  \
%                              |     \
%                              |      \
%                              |       \
%                              |        \
%                              | *     * \
%                              | 1     2  \
%                              +-----+-----+ ----> r
%
%% Class definition
%
classdef Gauss_Tria < fem.Gauss
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Gauss_Tria()
            this = this@fem.Gauss(fem.Gauss.TRIANGULAR);
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
        function [ngp,w,gp] = quadrature(~,order)  
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
    end
end