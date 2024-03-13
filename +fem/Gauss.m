%% Gauss Class
%
%% Description
%
% This is an abstract super-class that generically specifies an
% <element.html element> Gauss integration point in the FEMOOLab program.
% Essentially, this super-class declares abstract methods that define the
% general behavior of Gauss quadrature for numerical integration of element
% matrices and forcing vectors. These abstract methods are the functions
% that should be implemented in a derived sub-class that deals with
% specific types of Gauss integration points.
%
%% Subclasses
%
% * <gauss_quad.html Gauss_Quad: quadrilateral quadrature subclass>
% * <gauss_tria.html Gauss_Tria: triangular quadrature subclass>
%
%% Class definition
%
classdef Gauss < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of quadrature
        QUADRILATERAL = int32(1); 
        TRIANGULAR    = int32(2);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;  % flag for type of Gauss quadrature
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
        % given order of line quadrature.
        % Output:
        %  ngp: number of Gauss points
        %  w:   vector of Gauss points weights
        %  gp:  array of Gauss points parametric coordinates
        function [ngp,w,gp] = lineQuadrature(~,order)
            % Number of Gauss points
            ngp = order;
            
            % Assemble arrays of weights and parametric coordinates
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