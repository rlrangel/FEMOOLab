%% Gauss_Quad Class (Quadrilateral Gauss Quadrature)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <gauss.html Gauss: integration quadrature super-class>
% to deal with quadrilateral element shapes.
%
% Eg.: Quadrature of order 2:
%                                  ^ s
%                                  |
%                                  |
%                           +---------------+
%                           |   3   |   4   |
%                           |   *   |   *   |
%                           |       |       |
%                           |       ------- | --> r
%                           |               |
%                           |   *       *   |
%                           |   1       2   |
%                           +---------------+
%
%% Class definition
%
classdef Gauss_Quad < fem.Gauss
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Gauss_Quad()
            this = this@fem.Gauss(fem.Gauss.QUADRILATERAL);
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
    end
end