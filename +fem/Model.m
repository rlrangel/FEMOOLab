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
        anm             = [];               % object of Anm (Analysis Model) class
        anme fem.Anme   = fem.Anme.empty    % object of Anme (Analysis Method) class
        res  drv.Result = drv.Result.empty; % object of Result class
        
        % Model properties
        nnp       int32        = int32.empty;        % number of nodes
        nodes                  = [];                 % vector of objects of Node class or CtrlPt class
        nel       int32        = int32.empty;        % number of elements
        elems                  = [];                 % vector of objects of Element_Isoparametric class
        nmat      int32        = int32.empty;        % number of materials
        materials fem.Material = fem.Material.empty; % vector of objects of Material class
        nep       int32        = int32.empty         % number of extrapolation nodes
        extNodes  fem.ExtNode  = fem.ExtNode.empty   % vector of objects of Node class (extrapolation nodes). 
                                                     % In isoparametric analysis, nodes and extrapolation nodes 
                                                     % are the same. 
        
        % Model properties for isogeometric analysis
        nsurf     int32        = int32.empty;        % number of surfaces
        surfaces  fem.Surface  = fem.Surface.empty   % vector of objects of Surface class
        
        % Degree-of-freedom numbering
        ID   int32 = int32.empty; % global d.o.f. numbering matrix
        neq  int32 = int32.empty; % number of d.o.f.'s
        neqf int32 = int32.empty; % number of free d.o.f.'s
        neqc int32 = int32.empty; % number of constrained (fixed) d.o.f.'s
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Model()
            this.res = drv.Result();
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
        function setupDOFNum(this)
            % Dimension global d.o.f. numbering matrix
            this.ID = zeros(this.anm.ndof,this.nnp);
            
            % Initialize number of fixed d.o.f.'s
            this.neqc = 0;
            
            % Count number of fixed d.o.f.'s and setup ID matrix
            for i = 1:this.nnp
                for j = 1:this.anm.ndof
                    % Get d.o.f number
                    dof = this.anm.gla(j);
                    
                    % Check for fixed d.o.f
                    if (this.nodes(i).fixdDOF(dof))
                        this.neqc = this.neqc + 1;
                        this.ID(j,i) = 1;
                    end
                end
            end
            
            % Compute total number of free d.o.f.
            this.neqf = this.neq - this.neqc;
        end
        
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
        % Assemble global initial conditions matrix.
        % Set initial values only to free d.o.f's.
        %  IC(i,j) = initial value of time derivative j-1 of d.o.f. i.
        function IC = gblInitCondMtx(this)
            % Initialize initial conditions matrix
            IC = zeros(this.neqf,3);
            
            for i = 1:this.nnp
                for j = 1:this.anm.ndof
                    % Get d.o.f numbers
                    id  = this.ID(j,i);
                    dof = this.anm.gla(j);
                    
                    % Apply initial conditions to free d.o.f.'s
                    if (id <= this.neqf)
                        if (~isempty(this.nodes(i).initDOF))
                            IC(id,1) = this.nodes(i).initDOF(dof);
                        end
                        if (~isempty(this.nodes(i).initDOFt))
                            IC(id,2) = this.nodes(i).initDOFt(dof);
                        end
                        if (~isempty(this.nodes(i).initDOFtt))
                            IC(id,3) = this.nodes(i).initDOFtt(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add values of natural boundary conditions prescribed at nodal
        % points to global forcing vector.
        function F = addPointForce(this,F)
            for i = 1:this.nnp
                if (~isempty(this.nodes(i).prscNBC))
                    for j = 1:this.anm.ndof
                        % Get d.o.f numbers
                        id  = this.ID(j,i);
                        dof = this.anm.gla(j);
                        
                        % Add prescribed NBC value to global vector
                        F(id) = F(id) + this.nodes(i).prscNBC(dof);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add equivalent nodal forcing contributions from natural boundary
        % conditions and internal source to global forcing vector.
        function F = addEquivForce(this,F)
            for i = 1:this.nel
                gle = this.elems(i).gle;
                
                % Add equivalent forcing vector from NBC prescribed over edges
                if (~isempty(this.elems(i).lineNBC1))
                    Fedge = this.elems(i).edgeEquivForceVct(this);
                    F(gle) = F(gle) + Fedge;
                end
                
                % Add equivalent forcing vector from internal domain source
                if (~isempty(this.elems(i).src))
                    Fdom = this.elems(i).domainEquivForceVct();
                    F(gle) = F(gle) + Fdom;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add essencial boundary conditions (prescribed d.o.f.'s values)
        % to global vector of state variables.
        % Set prescribed values only to fixed d.o.f's.
        function U = addPrescDOF(this,U)
            for i = 1:this.nnp
                if (~isempty(this.nodes(i).prscDOF))
                    for j = 1:this.anm.ndof
                        % Get d.o.f. numbers
                        id  = this.ID(j,i);
                        dof = this.anm.gla(j);
                        
                        % Add prescribed d.o.f. value to global vector
                        if (this.nodes(i).fixdDOF(dof))
                            U(id) = this.nodes(i).prscDOF(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute maximum number nodes of all elements.
        function nep = maxNumElemExtNodes(this)
            nep = 0;
            for i = 1:this.nel
                if this.elems(i).shape.nep > nep
                    nep = this.elems(i).shape.nep;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute maximum number of derived variable Gauss points of all elements.
        function npts = maxGaussDerivedVarNpts(this)
            npts = 0;
            for i = 1:this.nel
                if this.elems(i).gderive_npts > npts
                    npts = this.elems(i).gderive_npts;
                end
            end
        end
    end
end