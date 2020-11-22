%% Model Class
%
%% Description
%
% This class defines a model object in the StAnOOP program.
% A model object is responsible for storing the global variables of a
% structural model.
% The methods of a model object are general functions that are not dependent
% on the <anm.html Anm: analysis model> or <shape.html Shape: element shape>
% type.
%
%% Main properties
%
%%%
% * Object of <anm.html Anm: analysis model class>
% * Object <result.html Result: Plot results driver class>
% * Vector of objects of <node.html Node: node class>
% * Vector of objects of <element.html Element: finite element class>
% * Vector of objects of <material.html Material: material class>
%
classdef Model < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General:
        anm        = [];                              % object of Anm (analysis model) class
        res        drv.Result = drv.Result.empty;     % object of Result (plot results driver) class
        
        % Model properties:
        nnp        int32 = int32.empty;               % number of nodes
        nodes      fem.Node = fem.Node.empty;         % vector of objects of Node class
        nel        int32 = int32.empty;               % number of elements
        elems      fem.Element = fem.Element.empty;   % vector of objects of Element (finite element) class
        nmat       int32 = int32.empty;               % number of materials
        materials  fem.Material = fem.Material.empty; % vector of objects of Material class
        nthks      int32 = int32.empty;               % number of thickness values
        thickness  double = double.empty;             % vector of element thicknesses
        nintord    int32 = int32.empty;               % number of integration orders
        intgrorder int32 = int32.empty;               % vector of pairs of integration orders
        
        % Degree-of-freedom numbering:
        ID         int32 = int32.empty;               % global d.o.f. numbering matrix
        neq        int32 = int32.empty;               % number of d.o.f.'s
        neqf       int32 = int32.empty;               % number of free d.o.f.'s
        neqc       int32 = int32.empty;               % number of constrained (fixed) d.o.f.'s
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function mdl = Model()
            return;
        end
    end
    
    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute stresses at Gauss points and principal stresses of all
        % elements.
        function gaussStressInplane(mdl,anl,res)
            maxElmGaussPts = mdl.maxGaussStressNpts;
            
            res.ngp     = zeros(mdl.nel,1);
            res.x_gp    = zeros(maxElmGaussPts*mdl.nel,1);
            res.y_gp    = zeros(maxElmGaussPts*mdl.nel,1);
            res.sxx_gp  = zeros(maxElmGaussPts,mdl.nel);
            res.syy_gp  = zeros(maxElmGaussPts,mdl.nel);
            res.txy_gp  = zeros(maxElmGaussPts,mdl.nel);
            res.s1_gp   = zeros(maxElmGaussPts,mdl.nel);
            res.s2_gp   = zeros(maxElmGaussPts,mdl.nel);
            res.tmax_gp = zeros(maxElmGaussPts,mdl.nel);
            res.s1x_gp  = zeros(maxElmGaussPts*mdl.nel,1);
            res.s1y_gp  = zeros(maxElmGaussPts*mdl.nel,1);
            res.s2x_gp  = zeros(maxElmGaussPts*mdl.nel,1);
            res.s2y_gp  = zeros(maxElmGaussPts*mdl.nel,1);
            
            npts = 0;
            for i = 1:mdl.nel
                d = res.D(mdl.elems(i).gle);
                
                [ngp,str,gpc] = mdl.elems(i).gaussStress(d);
                res.ngp(i)          = ngp;
                res.sxx_gp(1:ngp,i) = str(1,1:ngp);
                res.syy_gp(1:ngp,i) = str(2,1:ngp);
                res.txy_gp(1:ngp,i) = str(3,1:ngp);
                
                for j = 1:ngp
                    npts = npts + 1;
                    res.x_gp(npts) = gpc(1,j);
                    res.y_gp(npts) = gpc(2,j);
                    
                    [prc,thetap] = anl.princStress(str(:,j));
                    res.s1_gp(j,i)   = prc(1);
                    res.s2_gp(j,i)   = prc(2);
                    res.tmax_gp(j,i) = prc(3);
                    res.s1x_gp(npts) = res.s1_gp(j,i)*cos(thetap);
                    res.s1y_gp(npts) = res.s1_gp(j,i)*sin(thetap);
                    res.s2x_gp(npts) = res.s2_gp(j,i)*cos(thetap+(pi/2.0));
                    res.s2y_gp(npts) = res.s2_gp(j,i)*sin(thetap+(pi/2.0));
                end
            end
            
            res.sxx_gp_min  = min(min(res.sxx_gp));
            res.sxx_gp_max  = max(max(res.sxx_gp));
            res.syy_gp_min  = min(min(res.syy_gp));
            res.syy_gp_max  = max(max(res.syy_gp));
            res.txy_gp_min  = min(min(res.txy_gp));
            res.txy_gp_max  = max(max(res.txy_gp));
            res.s1_gp_min   = min(min(res.s1_gp));
            res.s1_gp_max   = max(max(res.s1_gp));
            res.s2_gp_min   = min(min(res.s2_gp));
            res.s2_gp_max   = max(max(res.s2_gp));
            res.tmax_gp_min = min(min(res.tmax_gp));
            res.tmax_gp_max = max(max(res.tmax_gp));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute node extrapolated stress components and principal stresses
        % for all elements.
        % The nodal stress components are computed by extrapolation of Gauss
        % point stress components using the TGN matrix.
        function elemStressExtrapInplane(mdl,res)
            maxNen = mdl.maxNumElemNodes;
            
            res.sxx_elemextrap  = zeros(maxNen,mdl.nel);
            res.syy_elemextrap  = zeros(maxNen,mdl.nel);
            res.txy_elemextrap  = zeros(maxNen,mdl.nel);
            res.s1_elemextrap   = zeros(maxNen,mdl.nel);
            res.s2_elemextrap   = zeros(maxNen,mdl.nel);
            res.tmax_elemextrap = zeros(maxNen,mdl.nel);
            
            for i = 1:mdl.nel
                TGN = mdl.elems(i).TGN;
                nen = mdl.elems(i).shape.nen;
                ngp = res.ngp(i);
                res.sxx_elemextrap(1:nen,i)  = TGN * res.sxx_gp(1:ngp,i);
                res.syy_elemextrap(1:nen,i)  = TGN * res.syy_gp(1:ngp,i);
                res.txy_elemextrap(1:nen,i)  = TGN * res.txy_gp(1:ngp,i);
                res.s1_elemextrap(1:nen,i)   = TGN * res.s1_gp(1:ngp,i);
                res.s2_elemextrap(1:nen,i)   = TGN * res.s2_gp(1:ngp,i);
                res.tmax_elemextrap(1:nen,i) = TGN * res.tmax_gp(1:ngp,i);
            end
            
            res.sxx_elemextrap_min  = min(min(res.sxx_elemextrap));
            res.sxx_elemextrap_max  = max(max(res.sxx_elemextrap));
            res.syy_elemextrap_min  = min(min(res.syy_elemextrap));
            res.syy_elemextrap_max  = max(max(res.syy_elemextrap));
            res.txy_elemextrap_min  = min(min(res.txy_elemextrap));
            res.txy_elemextrap_max  = max(max(res.txy_elemextrap));
            res.s1_elemextrap_min   = min(min(res.s1_elemextrap));
            res.s1_elemextrap_max   = max(max(res.s1_elemextrap));
            res.s2_elemextrap_min   = min(min(res.s2_elemextrap));
            res.s2_elemextrap_max   = max(max(res.s2_elemextrap));
            res.tmax_elemextrap_min = min(min(res.tmax_elemextrap));
            res.tmax_elemextrap_max = max(max(res.tmax_elemextrap));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute extrapolated node smoothed stress components and principal
        % stresses for all nodes.
        % The nodal stress components are computed by averaging values of
        % element extrapolated nodal stress components of all elements
        % adjacent to each node.
        function nodeStressExtrapInplane(mdl,res)
            
            % assemble vector with number of adjacent elements of each node
            node_adjelems = zeros(mdl.nnp,1);
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    n = mdl.elems(i).shape.nodes(j).id;
                    node_adjelems(n) = node_adjelems(n) + 1;
                end
            end
            
            res.sxx_nodeextrap  = zeros(mdl.nnp,1);
            res.syy_nodeextrap  = zeros(mdl.nnp,1);
            res.txy_nodeextrap  = zeros(mdl.nnp,1);
            res.s1_nodeextrap   = zeros(mdl.nnp,1);
            res.s2_nodeextrap   = zeros(mdl.nnp,1);
            res.tmax_nodeextrap = zeros(mdl.nnp,1);
            
            for i = 1:mdl.nel
                nen = mdl.elems(i).shape.nen;
                for j = 1:nen
                    n = mdl.elems(i).shape.nodes(j).id;
                    res.sxx_nodeextrap(n)  = res.sxx_nodeextrap(n)  + res.sxx_elemextrap(j,i);
                    res.syy_nodeextrap(n)  = res.syy_nodeextrap(n)  + res.syy_elemextrap(j,i);
                    res.txy_nodeextrap(n)  = res.txy_nodeextrap(n)  + res.txy_elemextrap(j,i);
                    res.s1_nodeextrap(n)   = res.s1_nodeextrap(n)   + res.s1_elemextrap(j,i);
                    res.s2_nodeextrap(n)   = res.s2_nodeextrap(n)   + res.s2_elemextrap(j,i);
                    res.tmax_nodeextrap(n) = res.tmax_nodeextrap(n) + res.tmax_elemextrap(j,i);
                end
            end
            
            for i = 1:mdl.nnp
                res.sxx_nodeextrap(i)  = res.sxx_nodeextrap(i)  / node_adjelems(i);
                res.syy_nodeextrap(i)  = res.syy_nodeextrap(i)  / node_adjelems(i);
                res.txy_nodeextrap(i)  = res.txy_nodeextrap(i)  / node_adjelems(i);
                res.s1_nodeextrap(i)   = res.s1_nodeextrap(i)   / node_adjelems(i);
                res.s2_nodeextrap(i)   = res.s2_nodeextrap(i)   / node_adjelems(i);
                res.tmax_nodeextrap(i) = res.tmax_nodeextrap(i) / node_adjelems(i);
            end
            
            res.sxx_nodeextrap_min  = min(min(res.sxx_nodeextrap));
            res.sxx_nodeextrap_max  = max(max(res.sxx_nodeextrap));
            res.syy_nodeextrap_min  = min(min(res.syy_nodeextrap));
            res.syy_nodeextrap_max  = max(max(res.syy_nodeextrap));
            res.txy_nodeextrap_min  = min(min(res.txy_nodeextrap));
            res.txy_nodeextrap_max  = max(max(res.txy_nodeextrap));
            res.s1_nodeextrap_min   = min(min(res.s1_nodeextrap));
            res.s1_nodeextrap_max   = max(max(res.s1_nodeextrap));
            res.s2_nodeextrap_min   = min(min(res.s2_nodeextrap));
            res.s2_nodeextrap_max   = max(max(res.s2_nodeextrap));
            res.tmax_nodeextrap_min = min(min(res.tmax_nodeextrap));
            res.tmax_nodeextrap_max = max(max(res.tmax_nodeextrap));
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Assemble global d.o.f numbering matrix:
        %  ID(j,i) = equation number of d.o.f. j of node i
        %  Free d.o.f.'s have the initial numbering.
        %  countF --> counts free d.o.f.'s. (numbered first)
        %  countC --> counts fixed d.o.f.'s. (numbered later)
        function assembleDOFNum(mdl)
            % Initialize numbers for free and fixed d.o.f.'s.
            countF = 0;
            countC = mdl.neqf;
            
            % Check if each d.o.f. is free or fixed to increment equation
            % number and store it in the ID matrix.
            for i = 1:mdl.nnp
                for j = 1:mdl.anm.ndof
                    if (mdl.ID(j,i) == 0)
                        countF = countF + 1;
                        mdl.ID(j,i) = countF;
                    else
                        countC = countC + 1;
                        mdl.ID(j,i) = countC;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assemble element gather vector (gle) that stores element global
        % d.o.f. numbers.
        function assembleGle(mdl)
            for i = 1:mdl.nel
                % Initialize element gather vector
                mdl.elems(i).gle = zeros(mdl.elems(i).shape.nen*mdl.anm.ndof,1);
                
                % Assemble d.o.f.'s of each node to element gather vector
                m = 0;
                for j = 1:mdl.elems(i).shape.nen
                    for k = 1:mdl.anm.ndof
                        m = m + 1;
                        mdl.elems(i).gle(m) = mdl.ID(k,mdl.elems(i).shape.nodes(j).id);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global elastic stiffness matrix.
        function Ke = gblElastStiffMtx(mdl)
            % Initialize global elastic matrix
            Ke = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element elastic matrix
                ke = mdl.elems(i).elastStiffMtx();
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                Ke(gle,gle) = Ke(gle,gle) + ke;
            end
        end
        
        %------------------------------------------------------------------
        % Add equivalent nodal loads from distributed line and area loads
        % to global forcing vector.
        function F = addEquivLoad(mdl,F)
            for i = 1:mdl.nel
                if (~isempty(mdl.elems(i).lineLoad))
                    % Get element equivalent nodal load vectors
                    fline = mdl.elems(i).edgeEquivLoadVct();
                    
                    % Assemble element vector to global vector
                    gle = mdl.elems(i).gle;
                    F(gle) = F(gle) + fline;
                end
                if (~isempty(mdl.elems(i).domainLoad))
                    % Get element equivalent nodal load vectors
                    farea = mdl.elems(i).domainEquivLoadVct();
                    
                    % Assemble element vector to global vector
                    gle = mdl.elems(i).gle;
                    F(gle) = F(gle) + farea;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute maximum number nodes of all elements.
        function nen = maxNumElemNodes(mdl)
            nen = 0;
            for i = 1:mdl.nel
                if mdl.elems(i).shape.nen > nen
                    nen = mdl.elems(i).shape.nen;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute maximum number of stress Gauss points of all elements.
        function npts = maxGaussStressNpts(mdl)
            npts = 0;
            for i = 1:mdl.nel
                if mdl.elems(i).gstress_npts > npts
                    npts = mdl.elems(i).gstress_npts;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute stresses at Gauss points and principal stresses of all
        % elements.
        function gaussStress(mdl,anl,res)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                mdl.gaussStressInplane(anl,res);
            end
        end
        
        %------------------------------------------------------------------
        % Compute node extrapolated stress components and principal stresses
        % for all elements.
        % The nodal stress components are computed by extrapolation of Gauss
        % point stress components using the TGN matrix.
        function elemStressExtrap(mdl,res)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                mdl.elemStressExtrapInplane(res);
            end
        end
        
        %------------------------------------------------------------------
        % Compute extrapolated node smoothed stress components and principal
        % stresses for all nodes.
        % The nodal stress components are computed by averaging values of
        % element extrapolated nodal stress components of all elements
        % adjacent to each node.
        function nodeStressExtrap(mdl,res)
            if mdl.anm.type == fem.Anm.PLANE_STRESS || ...
               mdl.anm.type == fem.Anm.PLANE_STRAIN || ...
               mdl.anm.type == fem.Anm.AXISYMMETRIC
                mdl.nodeStressExtrapInplane(res);
            end
        end
    end
end