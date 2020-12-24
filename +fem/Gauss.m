%% Gauss Class
%
%% Description
%
% This is an abstract super-class that generically specifies an
% <element.html element> Gauss integration point in the FEMOOLab program.
%
% Essentially, this super-class declares abstract methods that define the
% particular behavior of a Gauss quadrature for integrating the stiffness
% matrix of a finite element or for integrating the equivalente nodal
% load of a finite load element. These abstract methods are the
% functions that should be implemented in a derived sub-class that deals
% with specific types of Gauss integration points.
%
%% Current subclasses:
%
%%%
% * <gauss_quad.html Gauss_Quad: quadrilateral quadrature subclass>
% * <gauss_tria.html Gauss_Tria: triangular quadrature subclass>
%
classdef Gauss < handle
    %% Constant values
    properties (Constant = true, Access = public)
        QUADRILATERAL = int32(1); 
        TRIANGULAR    = int32(2);
    end
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Gauss points
        type      int32 = int32.empty;   % flag for type of Gauss quadrature
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Gauss(type)
            this.type = type;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Get parametric coordinates and weights of Gauss points for a
        % given quadrature order.
        [ngp,w,gp] = quadrature(this,order);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Get parametric coordinates and weights of edge Gauss points for a
        % line quadrature and given order.
        % Input:
        %  order: quadrature order
        % Output:
        %  ngp: number of Gauss points
        %  w:   vector of Gauss points weights
        %  gp:  array of Gauss points parametric coordinates
        function [ngp,w,gp] = lineQuadrature(~,order)
            % Number of Gauss points
            ngp = order;
            
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
            for i = 1:order
                w(i)    = qwgt(i,order);
                gp(1,i) = qpts(i,order);
                gp(2,i) = 0.0;
            end
        end
    end
end