%% Gauss_Lin Class
%
%% Description
%
% This is a sub-class of the <Gauss.html Gauss> class for the
% implementation of *Linear Quadrature*.
%
classdef Gauss_Lin < Gauss
    %% Constructor method
    methods
        function this = Gauss_Lin()
            this = this@Gauss(Gauss.LINEAR);
        end
    end
    
    %% Public methods: implementation of super-class declarations
    methods
        %------------------------------------------------------------------
        function [ngp,w,gp] = Quadrature(~,order)
            ngp = order;
            
            if (order == 1)
                w  = 2.00000000000000;

                gp = 0.000000000000000;

            elseif (order == 2)
                w  = [  1.000000000000000  1.000000000000000 ];

                gp = [ -0.577350269189626  0.577350269189626 ];

            elseif (order == 3)
                w  = [  0.555555555555555  0.888888888888888  0.555555555555555 ];

                gp = [ -0.774596669241483  0.000000000000000  0.774596669241483 ];

            elseif (order == 4)
                w  = [  0.347854845137454  0.652145154862546  0.652145154862546  0.347854845137454 ];
                
                gp = [ -0.861136311594053 -0.339981043584856  0.339981043584856  0.861136311594053 ];
            end
        end
    end
end