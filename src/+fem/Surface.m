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
        id            int32  = int32.empty;  % surface identification number
        ctrlNet       double = double.empty; % vector of id's of control points in surface
        elems         double = double.empty; % vector of id's of elements in surface
        knotVectorXi  double = double.empty; % knotvector in xi direction
        knotVectorEta double = double.empty; % knotvector in eta direction
        degreeXi      double = double.empty; % degree in xi direction
        degreeEta     double = double.empty; % degree in eta direction
        weights       double = double.empty; % weights of points in control net
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Surface()
            return;
        end
    end
end