%% Gauss_Tri Class
%
%% Description
%
% This is a sub-class of the <Gauss.html Gauss> class for the
% implementation of *Triangular Quadrature*.
%
% <<../images/tutorials/gauss_triangular.png>>
%
classdef Gauss_Tri < Gauss
    %% Constructor method
    methods
        function this = Gauss_Tri()
            this = this@Gauss(Gauss.TRIANGULAR);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function [ngp,w,gp] = Quadrature(~,order)
            if (order == 1)
                ngp = 1;

                w = 0.5;

                gp = [ 0.333333333333
                       0.333333333333 ];

            elseif (order == 2 || order == 3)
                ngp = 3;

                w = [ 0.166666666666 0.166666666666 0.166666666666 ];

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