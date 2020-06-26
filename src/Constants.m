%% Constants Class
%
% This class defines an object with constant flag variables as properties
% to be used in the StAnOOP program.
% An object of this class must be created in every function that uses a
% constant variable.
%
classdef Constants < handle
    %% Public properties (Access Only)
    properties (SetAccess = private, GetAccess = public)
        % Types of analysis
        LINEAR_ELASTIC = 1;
        
        % Types of analysis model
        PLANE_STRESS = 1;
        PLANE_STRAIN = 2;
        AXISYMMETRIC = 3;
        
        % Types of element
        TRI3  = 1;            % linear triangular element
        QUAD4 = 2;            % bilinear quadrilateral element
        TRI6  = 3;            % quadratic triangular element
        QUAD8 = 4;            % serendipity quadratic quadrilateral element
        
        % Types of finite element edge
        LINE2 = 1;            % linear edge        
        LINE3 = 2;            % quadratic edge
        
        % Types of Gauss quadrature
        QUADRATURE_LINE = 1;  % line quadrature
        QUADRATURE_TRIA = 2;  % triangular quadrature
        QUADRATURE_QUAD = 3;  % quadrilateral quadrature
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function c = Constants()
            return;
        end
    end
end