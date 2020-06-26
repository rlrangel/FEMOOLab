%% Model Class
%
% This class defines a model object in the StAnOOP program.
% A model object is responsible for storing the global variables of a
% structural model.
% The methods of a model object are general functions that are not dependent
% on the analysis model or element type.
%
classdef Model < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General:
        anm = [];         % object of the Anm (Analysis Model) class
        res = [];         % object of the Result class
        
        % Model properties:
        nnp       = 0;    % number of nodes
        nodes     = [];   % vector of objects of the Node class
        nel       = 0;    % number of elements
        elems     = [];   % vector of objects of the Element class
        nmat      = 0;    % number of materials
        materials = [];   % vector of objects of the Material class
        
        % Degree-of-freedom numbering:
        ID   = [];        % global d.o.f. numbering matrix
        neq  = 0;         % number of d.o.f.'s
        neqf = 0;         % number of free d.o.f.'s
        neqc = 0;         % number of constrained (fixed) d.o.f.'s
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function mdl = Model()
            return;
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Initialize global d.o.f. numbering matrix with ones and zeros,
        % and count total number of equations of free and fixed  d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(j,i) = 0, d.o.f. j of node i is free.
        %  if ID(j,i) = 1, d.o.f. j of node i is fixed.
        function setupDOFNum(mdl)
            % Dimension global d.o.f. numbering matrix
            mdl.ID = zeros(mdl.anm.ndof,mdl.nnp);
            
            % Initialize number of fixed d.o.f.'s
            mdl.neqc = 0;
            
            % Count number of fixed d.o.f.'s and setup ID matrix
            for i = 1:mdl.nnp
                for j = 1:mdl.anm.ndof
                    if (mdl.nodes(i).ebc(j) == 1)
                        mdl.neqc = mdl.neqc + 1;
                        mdl.ID(j,i) = 1;
                    end
                end
            end
            
            % Compute total number of free d.o.f.
            mdl.neqf = mdl.neq - mdl.neqc;
        end
        
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
                mdl.elems(i).gle = zeros(mdl.elems(i).nen*mdl.anm.ndof,1);
                
                % Assemble d.o.f.'s of each node to element gather vector
                m = 0;
                for j = 1:mdl.elems(i).nen
                    for k = 1:mdl.anm.ndof
                        m = m + 1;
                        mdl.elems(i).gle(m) = mdl.ID(k,mdl.elems(i).nodes(j).id);
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
        % Add point loads to global forcing vector, including the
        % components that correspond to fixed d.o.f.'s.
        function F = addPointLoad(mdl,F)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).load))
                    for j = 1:mdl.anm.ndof
                        id = mdl.ID(j,i);
                        F(id) = mdl.nodes(i).load(j);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add equivalent nodal loads from distributed line and area loads
        % to global forcing vector.
        function F = addEquivLoad(mdl,F)
            for i = 1:mdl.nel
                if (~isempty(mdl.elems(i).lineLoad))
                    % Get element equivalent nodal load vectors
                    fline = mdl.elems(i).lineEquivLoadVct();
                    
                    % Assemble element vector to global vector
                    gle = mdl.elems(i).gle;
                    F(gle) = F(gle) + fline;
                end
                if (~isempty(mdl.elems(i).areaLoad))
                    % Get element equivalent nodal load vectors
                    farea = mdl.elems(i).areaEquivLoadVct();
                    
                    % Assemble element vector to global vector
                    gle = mdl.elems(i).gle;
                    F(gle) = F(gle) + farea;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add nodal prescribed displacements (known support settlement values)
        % to global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of global displacement vector that corresponds to a free d.o.f.
        function D = addPrescDispl(mdl,D)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).prescDispl))
                    for j = 1:mdl.anm.ndof
                        id = mdl.ID(j,i);
                        if (id > mdl.neqf)
                            D(id) = mdl.nodes(i).prescDispl(j);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute stresses at Gauss points and principal stresses of all
        % elements.
        function gaussStress(mdl,anl,res)
            pts = mdl.elems(1).gstress_npts;  % TEMPORARY SIMPLIFICATION !!
            
            res.x_gp    = zeros(pts*mdl.nel,1);
            res.y_gp    = zeros(pts*mdl.nel,1);
            res.sx_gp   = zeros(pts,mdl.nel);
            res.sy_gp   = zeros(pts,mdl.nel);
            res.txy_gp  = zeros(pts,mdl.nel);
            res.s1_gp   = zeros(pts,mdl.nel);
            res.s2_gp   = zeros(pts,mdl.nel);
            res.tmax_gp = zeros(pts,mdl.nel);
            res.s1x_gp  = zeros(pts*mdl.nel,1);
            res.s1y_gp  = zeros(pts*mdl.nel,1);
            res.s2x_gp  = zeros(pts*mdl.nel,1);
            res.s2y_gp  = zeros(pts*mdl.nel,1);
            
            npts = 0;
            for i = 1:mdl.nel
                d = res.D(mdl.elems(i).gle);
                
                [ngp,str,gpc] = mdl.elems(i).gaussStress(d);
                res.sx_gp(:,i)  = str(1,:);
                res.sy_gp(:,i)  = str(2,:);
                res.txy_gp(:,i) = str(3,:);
                
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
            
            res.sx_gp_min   = min(min(res.sx_gp));
            res.sx_gp_max   = max(max(res.sx_gp));
            res.sy_gp_min   = min(min(res.sy_gp));
            res.sy_gp_max   = max(max(res.sy_gp));
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
        % Compute node extrapolated stress components and principal stresses
        % for all elements.
        % The nodal stress components are computed by extrapolation of Gauss
        % point stress components using the TGN matrix.
        function elemStressExtrap(mdl,res)
            nen = mdl.elems(1).nen;           % TEMPORARY SIMPLIFICATION !!
            
            res.sx_elemextrap   = zeros(nen,mdl.nel);
            res.sy_elemextrap   = zeros(nen,mdl.nel);
            res.txy_elemextrap  = zeros(nen,mdl.nel);
            res.s1_elemextrap   = zeros(nen,mdl.nel);
            res.s2_elemextrap   = zeros(nen,mdl.nel);
            res.tmax_elemextrap = zeros(nen,mdl.nel);
            
            for i = 1:mdl.nel
                TGN = mdl.elems(i).TGN;
                res.sx_elemextrap(:,i)   = TGN * res.sx_gp(:,i);
                res.sy_elemextrap(:,i)   = TGN * res.sy_gp(:,i);
                res.txy_elemextrap(:,i)  = TGN * res.txy_gp(:,i);
                res.s1_elemextrap(:,i)   = TGN * res.s1_gp(:,i);
                res.s2_elemextrap(:,i)   = TGN * res.s2_gp(:,i);
                res.tmax_elemextrap(:,i) = TGN * res.tmax_gp(:,i);
            end
            
            res.sx_elemextrap_min   = min(min(res.sx_elemextrap));
            res.sx_elemextrap_max   = max(max(res.sx_elemextrap));
            res.sy_elemextrap_min   = min(min(res.sy_elemextrap));
            res.sy_elemextrap_max   = max(max(res.sy_elemextrap));
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
        % Compute extrapolated node smoothed stress components and principal
        % stresses for all nodes.
        % The nodal stress components are computed by averaging values of
        % element extrapolated nodal stress components of all elements
        % adjacent to each node.
        function nodeStressExtrap(mdl,res)
            nen = mdl.elems(1).nen;           % TEMPORARY SIMPLIFICATION !!
            
            % assemble vector with number of adjacent elements of each node
            node_adjelems = zeros(mdl.nnp,1);
            for i = 1:mdl.nel
                for j = 1:nen
                    n = mdl.elems(i).nodes(j).id;
                    node_adjelems(n) = node_adjelems(n) + 1;
                end
            end
            
            res.sx_nodeextrap   = zeros(mdl.nnp,1);
            res.sy_nodeextrap   = zeros(mdl.nnp,1);
            res.txy_nodeextrap  = zeros(mdl.nnp,1);
            res.s1_nodeextrap   = zeros(mdl.nnp,1);
            res.s2_nodeextrap   = zeros(mdl.nnp,1);
            res.tmax_nodeextrap = zeros(mdl.nnp,1);
            
            for i = 1:mdl.nel
                for j = 1:nen
                    n = mdl.elems(i).nodes(j).id;
                    res.sx_nodeextrap(n)   = res.sx_nodeextrap(n)   + res.sx_elemextrap(j,i);
                    res.sy_nodeextrap(n)   = res.sy_nodeextrap(n)   + res.sy_elemextrap(j,i);
                    res.txy_nodeextrap(n)  = res.txy_nodeextrap(n)  + res.txy_elemextrap(j,i);
                    res.s1_nodeextrap(n)   = res.s1_nodeextrap(n)   + res.s1_elemextrap(j,i);
                    res.s2_nodeextrap(n)   = res.s2_nodeextrap(n)   + res.s2_elemextrap(j,i);
                    res.tmax_nodeextrap(n) = res.tmax_nodeextrap(n) + res.tmax_elemextrap(j,i);
                end
            end
            
            for i = 1:mdl.nnp
                res.sx_nodeextrap(i)   = res.sx_nodeextrap(i)   / node_adjelems(i);
                res.sy_nodeextrap(i)   = res.sy_nodeextrap(i)   / node_adjelems(i);
                res.txy_nodeextrap(i)  = res.txy_nodeextrap(i)  / node_adjelems(i);
                res.s1_nodeextrap(i)   = res.s1_nodeextrap(i)   / node_adjelems(i);
                res.s2_nodeextrap(i)   = res.s2_nodeextrap(i)   / node_adjelems(i);
                res.tmax_nodeextrap(i) = res.tmax_nodeextrap(i) / node_adjelems(i);
            end
            
            res.sx_nodeextrap_min   = min(min(res.sx_nodeextrap));
            res.sx_nodeextrap_max   = max(max(res.sx_nodeextrap));
            res.sy_nodeextrap_min   = min(min(res.sy_nodeextrap));
            res.sy_nodeextrap_max   = max(max(res.sy_nodeextrap));
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
end