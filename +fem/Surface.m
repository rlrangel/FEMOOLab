%% Surface Class
%
%% Description
%
%
%% Class definition
%
classdef Surface < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General properties
        id          int32  = int32.empty;  % surface identification number
        ctrlNet     double = double.empty; % vector of id's of control points in surface
        elems       double = double.empty; % vector of id's of elements in surface
        knotVectorU double = double.empty; % vector of knotvector in u direction
        knotVectorV double = double.empty; % vector of knotvector in v direction
        degreeU     double = double.empty; % degree in u direction
        degreeV     double = double.empty; % degree in v direction
        weights     double = double.empty; % weights of points in control net
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Surface()
            return;
        end
    end
end