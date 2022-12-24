%% Gauss Class
%
%% Description
%
% This is a handle super-class for the definition of an element Gauss
% quadrature for numerical integration of element arrays.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <Gauss_Lin.html Gauss_Lin>
% * <Gauss_Tri.html Gauss_Tri>
% * <Gauss_Quad.html Gauss_Quad>
%
classdef Gauss < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of quadrature
        LINEAR        = int32(1);
        TRIANGULAR    = int32(2);
        QUADRILATERAL = int32(3); 
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;   % flag for type of Gauss quadrature
    end
    
    %% Constructor method
    methods
        function this = Gauss(type)
            if (nargin > 0)
                this.type = type;
            end
        end
    end
    
    %% Abstract methods: implemented in derived sub-classes
    methods (Abstract)
        %------------------------------------------------------------------
        % Get parametric coordinates and weights of Gauss points for a
        % given quadrature order.
        % Output:
        %  ngp: number of Gauss points
        %  w:   vector of Gauss points weights
        %  gp:  array of Gauss points parametric coordinates
        [ngp,w,gp] = Quadrature(this,order);
    end
end