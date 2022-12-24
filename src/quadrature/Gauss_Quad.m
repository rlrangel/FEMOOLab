%% Gauss_Quad Class 
%
%% Description
%
% This is a sub-class of the <Gauss.html Gauss> class for the
% implementation of *Quadrilateral Quadrature*.
%
% <<../images/tutorials/gauss_quadrilateral.png>>
%
classdef Gauss_Quad < Gauss
    %% Constructor method
    methods
        function this = Gauss_Quad()
            this = this@Gauss(Gauss.QUADRILATERAL);
        end
    end
    
    %% Public methods: implementation of super-class declarations
    methods
        %------------------------------------------------------------------
        function [ngp,w,gp] = Quadrature(~,order)
            ngp = order * order;
            
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
            npts = 0;
            w    = zeros(1,ngp);
            gp   = zeros(2,ngp);
            for i = 1:order
                for j = 1:order
                    npts = npts + 1;
                    w(npts)    = qwgt(j,order) * qwgt(i,order);
                    gp(1,npts) = qpts(i,order);
                    gp(2,npts) = qpts(j,order);
                end
            end
        end
    end
end