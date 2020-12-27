%% Model Class
%
%% Description
%
% This class defines a model object in the FEMOOLab program.
% A model object is responsible for storing the global variables of a
% computational model.
% The methods of a model object are general functions that are not dependent
% on the analysis type, analysis model, or element shape.
%
%% Class definition
%
classdef Model < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        anm = [];                                     % object of Anm (Analysis Model) class
        res drv.Result = drv.Result.empty;            % object of Result class
        
        % Model properties
        nnp        int32        = int32.empty;        % number of nodes
        nodes      fem.Node     = fem.Node.empty;     % vector of objects of Node class
        nel        int32        = int32.empty;        % number of elements
        elems      fem.Element  = fem.Element.empty;  % vector of objects of Element (finite element) class
        nmat       int32        = int32.empty;        % number of materials
        materials  fem.Material = fem.Material.empty; % vector of objects of Material class
        
        % Degree-of-freedom numbering
        ID   int32 = int32.empty;                     % global d.o.f. numbering matrix
        neq  int32 = int32.empty;                     % number of d.o.f.'s
        neqf int32 = int32.empty;                     % number of free d.o.f.'s
        neqc int32 = int32.empty;                     % number of constrained (fixed) d.o.f.'s
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Model()
            this.res = drv.Result();
        end
    end
    
    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute stresses at Gauss points and principal stresses of all
        % elements.
        function gaussStressInplane(this,anl)
            r = this.res;
            maxElmGaussPts = this.maxGaussStressNpts;
            
            r.ngp     = zeros(this.nel,1);
            r.x_gp    = zeros(maxElmGaussPts*this.nel,1);
            r.y_gp    = zeros(maxElmGaussPts*this.nel,1);
            r.sxx_gp  = zeros(maxElmGaussPts,this.nel);
            r.syy_gp  = zeros(maxElmGaussPts,this.nel);
            r.txy_gp  = zeros(maxElmGaussPts,this.nel);
            r.s1_gp   = zeros(maxElmGaussPts,this.nel);
            r.s2_gp   = zeros(maxElmGaussPts,this.nel);
            r.tmax_gp = zeros(maxElmGaussPts,this.nel);
            r.s1x_gp  = zeros(maxElmGaussPts*this.nel,1);
            r.s1y_gp  = zeros(maxElmGaussPts*this.nel,1);
            r.s2x_gp  = zeros(maxElmGaussPts*this.nel,1);
            r.s2y_gp  = zeros(maxElmGaussPts*this.nel,1);
            
            npts = 0;
            for i = 1:this.nel
                u = r.U(this.elems(i).gle);
                
                [ngp,str,gpc] = this.elems(i).gaussStress(u);
                r.ngp(i)          = ngp;
                r.sxx_gp(1:ngp,i) = str(1,1:ngp);
                r.syy_gp(1:ngp,i) = str(2,1:ngp);
                r.txy_gp(1:ngp,i) = str(3,1:ngp);
                
                for j = 1:ngp
                    npts = npts + 1;
                    r.x_gp(npts) = gpc(1,j);
                    r.y_gp(npts) = gpc(2,j);
                    
                    % Compute principal stresses
                    [prc,thetap] = anl.princStress(str(:,j));
                    r.s1_gp(j,i)   = prc(1);
                    r.s2_gp(j,i)   = prc(2);
                    r.tmax_gp(j,i) = prc(3);
                    r.s1x_gp(npts) = r.s1_gp(j,i)*cos(thetap);
                    r.s1y_gp(npts) = r.s1_gp(j,i)*sin(thetap);
                    r.s2x_gp(npts) = r.s2_gp(j,i)*cos(thetap+(pi/2.0));
                    r.s2y_gp(npts) = r.s2_gp(j,i)*sin(thetap+(pi/2.0));
                end
            end
            
            r.sxx_gp_min  = min(min(r.sxx_gp));
            r.sxx_gp_max  = max(max(r.sxx_gp));
            r.syy_gp_min  = min(min(r.syy_gp));
            r.syy_gp_max  = max(max(r.syy_gp));
            r.txy_gp_min  = min(min(r.txy_gp));
            r.txy_gp_max  = max(max(r.txy_gp));
            r.s1_gp_min   = min(min(r.s1_gp));
            r.s1_gp_max   = max(max(r.s1_gp));
            r.s2_gp_min   = min(min(r.s2_gp));
            r.s2_gp_max   = max(max(r.s2_gp));
            r.tmax_gp_min = min(min(r.tmax_gp));
            r.tmax_gp_max = max(max(r.tmax_gp));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute stresses at Gauss points and principal stresses of all
        % elements.
        function gaussFluxInplane(this,~)
            r = this.res;
            maxElmGaussPts = this.maxGaussStressNpts;
            
            r.ngp      = zeros(this.nel,1);
            r.x_gp     = zeros(maxElmGaussPts*this.nel,1);
            r.y_gp     = zeros(maxElmGaussPts*this.nel,1);
            r.fxx_gp   = zeros(maxElmGaussPts,this.nel);
            r.fyy_gp   = zeros(maxElmGaussPts,this.nel);
            r.fp_gp    = zeros(maxElmGaussPts,this.nel);
            r.fpx_gp   = zeros(maxElmGaussPts*this.nel,1);
            r.fpy_gp   = zeros(maxElmGaussPts*this.nel,1);
            
            npts = 0;
            for i = 1:this.nel
                u = r.U(this.elems(i).gle);
                
                [ngp,flx,gpc] = this.elems(i).gaussStress(u);
                r.ngp(i)          = ngp;
                r.fxx_gp(1:ngp,i) = flx(1,1:ngp);
                r.fyy_gp(1:ngp,i) = flx(2,1:ngp);
                
                for j = 1:ngp
                    npts = npts + 1;
                    r.x_gp(npts) = gpc(1,j);
                    r.y_gp(npts) = gpc(2,j);
                    
                    % Compute principal stresses
                    prc = sqrt(flx(1,j)^2+flx(2,j)^2);
                    thetap = atan2(flx(2,j),flx(1,j));
                    
                    r.fp_gp(j,i)   = prc;
                    r.fpx_gp(npts) = r.fp_gp(j,i)*cos(thetap);
                    r.fpy_gp(npts) = r.fp_gp(j,i)*sin(thetap);
                end
            end
            
            r.fxx_gp_min = min(min(r.fxx_gp));
            r.fxx_gp_max = max(max(r.fxx_gp));
            r.fyy_gp_min = min(min(r.fyy_gp));
            r.fyy_gp_max = max(max(r.fyy_gp));
            r.fp_gp_min  = min(min(r.fp_gp));
            r.fp_gp_max  = max(max(r.fp_gp));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute node extrapolated stress components and principal stresses
        % for all elements.
        % The nodal stress components are computed by extrapolation of Gauss
        % point stress components using the TGN matrix.
        function elemStressExtrapInplane(this)
            r = this.res;
            maxNen = this.maxNumElemNodes;
            
            r.sxx_elemextrap  = zeros(maxNen,this.nel);
            r.syy_elemextrap  = zeros(maxNen,this.nel);
            r.txy_elemextrap  = zeros(maxNen,this.nel);
            r.s1_elemextrap   = zeros(maxNen,this.nel);
            r.s2_elemextrap   = zeros(maxNen,this.nel);
            r.tmax_elemextrap = zeros(maxNen,this.nel);
            
            for i = 1:this.nel
                TGN = this.elems(i).TGN;
                nen = this.elems(i).shape.nen;
                ngp = r.ngp(i);
                r.sxx_elemextrap(1:nen,i)  = TGN * r.sxx_gp(1:ngp,i);
                r.syy_elemextrap(1:nen,i)  = TGN * r.syy_gp(1:ngp,i);
                r.txy_elemextrap(1:nen,i)  = TGN * r.txy_gp(1:ngp,i);
                r.s1_elemextrap(1:nen,i)   = TGN * r.s1_gp(1:ngp,i);
                r.s2_elemextrap(1:nen,i)   = TGN * r.s2_gp(1:ngp,i);
                r.tmax_elemextrap(1:nen,i) = TGN * r.tmax_gp(1:ngp,i);
            end
            
            r.sxx_elemextrap_min  = min(min(r.sxx_elemextrap));
            r.sxx_elemextrap_max  = max(max(r.sxx_elemextrap));
            r.syy_elemextrap_min  = min(min(r.syy_elemextrap));
            r.syy_elemextrap_max  = max(max(r.syy_elemextrap));
            r.txy_elemextrap_min  = min(min(r.txy_elemextrap));
            r.txy_elemextrap_max  = max(max(r.txy_elemextrap));
            r.s1_elemextrap_min   = min(min(r.s1_elemextrap));
            r.s1_elemextrap_max   = max(max(r.s1_elemextrap));
            r.s2_elemextrap_min   = min(min(r.s2_elemextrap));
            r.s2_elemextrap_max   = max(max(r.s2_elemextrap));
            r.tmax_elemextrap_min = min(min(r.tmax_elemextrap));
            r.tmax_elemextrap_max = max(max(r.tmax_elemextrap));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute node extrapolated stress components and principal stresses
        % for all elements.
        % The nodal stress components are computed by extrapolation of Gauss
        % point stress components using the TGN matrix.
        function elemFluxExtrapInplane(this)
            r = this.res;
            maxNen = this.maxNumElemNodes;
            
            r.fxx_elemextrap = zeros(maxNen,this.nel);
            r.fyy_elemextrap = zeros(maxNen,this.nel);
            r.fp_elemextrap  = zeros(maxNen,this.nel);
            
            for i = 1:this.nel
                TGN = this.elems(i).TGN;
                nen = this.elems(i).shape.nen;
                ngp = r.ngp(i);
                r.fxx_elemextrap(1:nen,i) = TGN * r.fxx_gp(1:ngp,i);
                r.fyy_elemextrap(1:nen,i) = TGN * r.fyy_gp(1:ngp,i);
                r.fp_elemextrap(1:nen,i)  = TGN * r.fp_gp(1:ngp,i);
            end
            
            r.fxx_elemextrap_min = min(min(r.fxx_elemextrap));
            r.fxx_elemextrap_max = max(max(r.fxx_elemextrap));
            r.fyy_elemextrap_min = min(min(r.fyy_elemextrap));
            r.fyy_elemextrap_max = max(max(r.fyy_elemextrap));
            r.fp_elemextrap_min  = min(min(r.fp_elemextrap));
            r.fp_elemextrap_max  = max(max(r.fp_elemextrap));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute extrapolated node smoothed stress components and principal
        % stresses for all nodes.
        % The nodal stress components are computed by averaging values of
        % element extrapolated nodal stress components of all elements
        % adjacent to each node.
        function nodeStressExtrapInplane(this)
            r = this.res;
            
            % assemble vector with number of adjacent elements of each node
            node_adjelems = zeros(this.nnp,1);
            for i = 1:this.nel
                nen = this.elems(i).shape.nen;
                for j = 1:nen
                    n = this.elems(i).shape.nodes(j).id;
                    node_adjelems(n) = node_adjelems(n) + 1;
                end
            end
            
            r.sxx_nodeextrap  = zeros(this.nnp,1);
            r.syy_nodeextrap  = zeros(this.nnp,1);
            r.txy_nodeextrap  = zeros(this.nnp,1);
            r.s1_nodeextrap   = zeros(this.nnp,1);
            r.s2_nodeextrap   = zeros(this.nnp,1);
            r.tmax_nodeextrap = zeros(this.nnp,1);
            
            for i = 1:this.nel
                nen = this.elems(i).shape.nen;
                for j = 1:nen
                    n = this.elems(i).shape.nodes(j).id;
                    r.sxx_nodeextrap(n)  = r.sxx_nodeextrap(n)  + r.sxx_elemextrap(j,i);
                    r.syy_nodeextrap(n)  = r.syy_nodeextrap(n)  + r.syy_elemextrap(j,i);
                    r.txy_nodeextrap(n)  = r.txy_nodeextrap(n)  + r.txy_elemextrap(j,i);
                    r.s1_nodeextrap(n)   = r.s1_nodeextrap(n)   + r.s1_elemextrap(j,i);
                    r.s2_nodeextrap(n)   = r.s2_nodeextrap(n)   + r.s2_elemextrap(j,i);
                    r.tmax_nodeextrap(n) = r.tmax_nodeextrap(n) + r.tmax_elemextrap(j,i);
                end
            end
            
            for i = 1:this.nnp
                r.sxx_nodeextrap(i)  = r.sxx_nodeextrap(i)  / node_adjelems(i);
                r.syy_nodeextrap(i)  = r.syy_nodeextrap(i)  / node_adjelems(i);
                r.txy_nodeextrap(i)  = r.txy_nodeextrap(i)  / node_adjelems(i);
                r.s1_nodeextrap(i)   = r.s1_nodeextrap(i)   / node_adjelems(i);
                r.s2_nodeextrap(i)   = r.s2_nodeextrap(i)   / node_adjelems(i);
                r.tmax_nodeextrap(i) = r.tmax_nodeextrap(i) / node_adjelems(i);
            end
            
            r.sxx_nodeextrap_min  = min(min(r.sxx_nodeextrap));
            r.sxx_nodeextrap_max  = max(max(r.sxx_nodeextrap));
            r.syy_nodeextrap_min  = min(min(r.syy_nodeextrap));
            r.syy_nodeextrap_max  = max(max(r.syy_nodeextrap));
            r.txy_nodeextrap_min  = min(min(r.txy_nodeextrap));
            r.txy_nodeextrap_max  = max(max(r.txy_nodeextrap));
            r.s1_nodeextrap_min   = min(min(r.s1_nodeextrap));
            r.s1_nodeextrap_max   = max(max(r.s1_nodeextrap));
            r.s2_nodeextrap_min   = min(min(r.s2_nodeextrap));
            r.s2_nodeextrap_max   = max(max(r.s2_nodeextrap));
            r.tmax_nodeextrap_min = min(min(r.tmax_nodeextrap));
            r.tmax_nodeextrap_max = max(max(r.tmax_nodeextrap));
        end
        
        %------------------------------------------------------------------
        % 2D inplane analysis model:
        % Compute extrapolated node smoothed stress components and principal
        % stresses for all nodes.
        % The nodal stress components are computed by averaging values of
        % element extrapolated nodal stress components of all elements
        % adjacent to each node.
        function nodeFluxExtrapInplane(this)
            r = this.res;
            
            % assemble vector with number of adjacent elements of each node
            node_adjelems = zeros(this.nnp,1);
            for i = 1:this.nel
                nen = this.elems(i).shape.nen;
                for j = 1:nen
                    n = this.elems(i).shape.nodes(j).id;
                    node_adjelems(n) = node_adjelems(n) + 1;
                end
            end
            
            r.fxx_nodeextrap = zeros(this.nnp,1);
            r.fyy_nodeextrap = zeros(this.nnp,1);
            r.fp_nodeextrap  = zeros(this.nnp,1);
            
            for i = 1:this.nel
                nen = this.elems(i).shape.nen;
                for j = 1:nen
                    n = this.elems(i).shape.nodes(j).id;
                    r.fxx_nodeextrap(n) = r.fxx_nodeextrap(n) + r.fxx_elemextrap(j,i);
                    r.fyy_nodeextrap(n) = r.fyy_nodeextrap(n) + r.fyy_elemextrap(j,i);
                    r.fp_nodeextrap(n)  = r.fp_nodeextrap(n)  + r.fp_elemextrap(j,i);
                end
            end
            
            for i = 1:this.nnp
                r.fxx_nodeextrap(i) = r.fxx_nodeextrap(i) / node_adjelems(i);
                r.fyy_nodeextrap(i) = r.fyy_nodeextrap(i) / node_adjelems(i);
                r.fp_nodeextrap(i)  = r.fp_nodeextrap(i)  / node_adjelems(i);
            end
            
            r.fxx_nodeextrap_min = min(min(r.fxx_nodeextrap));
            r.fxx_nodeextrap_max = max(max(r.fxx_nodeextrap));
            r.fyy_nodeextrap_min = min(min(r.fyy_nodeextrap));
            r.fyy_nodeextrap_max = max(max(r.fyy_nodeextrap));
            r.fp_nodeextrap_min  = min(min(r.fp_nodeextrap));
            r.fp_nodeextrap_max  = max(max(r.fp_nodeextrap));
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
        function assembleDOFNum(this)
            % Initialize numbers for free and fixed d.o.f.'s.
            countF = 0;
            countC = this.neqf;
            
            % Check if each d.o.f. is free or fixed to increment equation
            % number and store it in the ID matrix.
            for i = 1:this.nnp
                for j = 1:this.anm.ndof
                    if (this.ID(j,i) == 0)
                        countF = countF + 1;
                        this.ID(j,i) = countF;
                    else
                        countC = countC + 1;
                        this.ID(j,i) = countC;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assemble element gather vector (gle) that stores element global
        % d.o.f. numbers.
        function assembleGle(this)
            for i = 1:this.nel
                % Initialize element gather vector
                this.elems(i).gle = zeros(this.elems(i).shape.nen*this.anm.ndof,1);
                
                % Assemble d.o.f.'s of each node to element gather vector
                m = 0;
                for j = 1:this.elems(i).shape.nen
                    for k = 1:this.anm.ndof
                        m = m + 1;
                        this.elems(i).gle(m) = this.ID(k,this.elems(i).shape.nodes(j).id);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute maximum number nodes of all elements.
        function nen = maxNumElemNodes(this)
            nen = 0;
            for i = 1:this.nel
                if this.elems(i).shape.nen > nen
                    nen = this.elems(i).shape.nen;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute maximum number of stress Gauss points of all elements.
        function npts = maxGaussStressNpts(this)
            npts = 0;
            for i = 1:this.nel
                if this.elems(i).gstress_npts > npts
                    npts = this.elems(i).gstress_npts;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute stresses at Gauss points and principal stresses of all
        % elements.
        function gaussStress(this,anl)
            if this.anm.type == fem.Anm.PLANE_STRESS || ...
               this.anm.type == fem.Anm.PLANE_STRAIN || ...
               this.anm.type == fem.Anm.AXISYM_STRESS
                this.gaussStressInplane(anl);
            elseif this.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                   this.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.gaussFluxInplane(anl);
            end
        end
        
        %------------------------------------------------------------------
        % Compute node extrapolated stress components and principal stresses
        % for all elements.
        % The nodal stress components are computed by extrapolation of Gauss
        % point stress components using the TGN matrix.
        function elemStressExtrap(this)
            if this.anm.type == fem.Anm.PLANE_STRESS || ...
               this.anm.type == fem.Anm.PLANE_STRAIN || ...
               this.anm.type == fem.Anm.AXISYM_STRESS
                this.elemStressExtrapInplane();
            elseif this.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                   this.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.elemFluxExtrapInplane();
            end
        end
        
        %------------------------------------------------------------------
        % Compute extrapolated node smoothed stress components and principal
        % stresses for all nodes.
        % The nodal stress components are computed by averaging values of
        % element extrapolated nodal stress components of all elements
        % adjacent to each node.
        function nodeStressExtrap(this)
            if this.anm.type == fem.Anm.PLANE_STRESS || ...
               this.anm.type == fem.Anm.PLANE_STRAIN || ...
               this.anm.type == fem.Anm.AXISYM_STRESS
                this.nodeStressExtrapInplane();
            elseif this.anm.type == fem.Anm.PLANE_CONDUCTION || ...
                   this.anm.type == fem.Anm.AXISYM_CONDUCTION
                this.nodeFluxExtrapInplane();
            end
        end
    end
end