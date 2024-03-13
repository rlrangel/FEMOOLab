%% Model class
classdef Model < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        anm = [];         % handle to an object of the Anm class
        elems = [];       % vector of handles to objects of the Elem class
        nodes = [];       % vector of handles to objects of the Node class
        materials = [];   % vector of handles to objects of the Material class
        nel = 0;          % number of elements
        nnp = 0;          % number of nodes
        nmat = 0;         % number of materials
        neq = 0;          % number of equations
        neqfree = 0;      % number of equations of free d.o.f.
        neqfixed = 0;     % number of equations of fixed d.o.f.
        ID = [];          % global d.o.f. numbering matrix
        K = [];           % global stiffness matrix
        F = [];           % global system forcing vector
        D = [];           % global system displacement vector
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function model = model(anm,elems,nodes,mat,nel,nnp,nmat,neq,neqfr,neqfx,K,F,D,ID)
            if (nargin > 0)
                model.anm = anm;
                model.elems = elems;
                model.nodes = nodes;
                model.materials = mat;
                model.nel = nel;
                model.nnp = nnp;
                model.nmat = nmat;
                model.neq = neq;
                model.neqfree = neqfr;
                model.neqfixed = neqfx;
                model.ID = ID;
                model.K = K;
                model.F = F;
                model.D = D;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Dimensions and initializes global stiffness matrix, global forcing
        % vector and global displacement vector.
        function dimKFD(model)
            model.K = zeros(model.neq,model.neq);
            model.F = zeros(model.neq,1);
            model.D = zeros(model.neq,1);
        end
        
        %------------------------------------------------------------------
        % Assembles global d.o.f (degree of freedom) numbering ID matrix:
        %  ID(k,n) = equation number of d.o.f. 'k' of node 'n'
        %  Free d.o.f.'s have the initial numbering.
        %  countF --> counts free or natural B.C. (unknown) d.o.f.
        %             (numbered first)
        %  countS --> counts constrainted by Support B.C.
        %             (essential - known) d.o.f.
        function assembleDOFNum(model)
            countF = 0;
            countS = model.neqfree;
            
            for n = 1:model.nnp
                for k = 1:model.anm.ndof
                    if(model.ID(k,n) == 0)
                        countF = countF + 1;
                        model.ID(k,n) = countF;
                    else
                        countS = countS + 1;
                        model.ID(k,n) = countS;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element gather vector (gle) that stores element d.o.f.
        % equation numbers.
        function assembleGle(model)
            for e = 1:model.nel
                model.elems(e).gle = zeros(model.elems(e).nen*model.anm.ndof,1);
                i = 0;
                
                for n = 1:model.elems(e).nen
                    for k = 1:model.anm.ndof
                        i = i + 1;
                        node = model.elems(e).nodes(n).id;
                        model.elems(e).gle(i) = model.ID(k,node);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles global stiffness matrix by adding the contribution of
        % all elements.
        function gblStiffMtx(model)
            for e = 1:model.nel
                ke = model.elems(e).stiffMtx();
                model.assembleElemMtx(ke,e);
            end
        end
        
        %------------------------------------------------------------------
        % Assembles given element stiffness matrix to global stiffness matrix.
        % Input arguments:
        %  keg: element local stiffness matrix in global system
        %  e: element identification number
        function assembleElemMtx(model,ke,e)
            gle = model.elems(e).gle;
            model.K(gle,gle) = model.K(gle,gle) + ke;
        end
        
        %------------------------------------------------------------------
        % Adds edge (element side) load as equivalent nodal loads (ENL) to
        % global forcing vector.
        function edgeLoads(model)
            for e = 1:model.nel
                for l = 1:size(model.elems(e).edgeLoad,1)
                    fe = model.elems(e).edgeENL(l);
                    model.assembleENL(fe,e);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Adds area (element) load as equivalent nodal loads (ENL) to
        % global forcing vector.
        function areaLoads(model)
            for e = 1:model.nel
                if isempty(model.elems(e).areaLoad) == 0
                    fe = model.elems(e).areaENL();
                    model.assembleENL(fe,e);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles given equivalent nodal load vector to any term of
        % the global forcing vector, including the terms that correspond
        % to constrained d.o.f.
        % Input arguments:
        %  fe: equivalent nodal load vector
        %  e: element identification number
        function assembleENL(model,fe,e)
            gle = model.elems(e).gle;
            model.F(gle) = model.F(gle) + fe;
        end
        
        %------------------------------------------------------------------
        % Partitions and solves the system of equations.
        % Partitions stiffness matrix K, forcing vector F and unknown vector D:
        %  f --> free or natural B.C. (unknown) degree-of-freedom
        %        (numbered first)
        %  c --> constrainted by support (essential - known) B.C. 
        %        degree-of-freedown
        %
        % [ Kff Kfc ] * [ Df ] = [ Ff ]
        % [ Kcf Kcc ]   [ Dc ] = [ Fc ]
        %
        % Output:
        %  stability: flag for stable free-free global matrix (0 -> unstable)
        function stability = solveEqnSystem(model)
            % Assemble free-free global matrix
            Kff = model.K(1:model.neqfree,1:model.neqfree);
            
            % Check for stable free-free global matrix by checking its
            % reciprocal condition number (a very low reciprocal condition
            % number indicates that the matrix is badly conditioned and
            % may be singular)
            stability = 1;
            if (rcond(Kff) < 10e-12)
                stability = 0;
                return;
            end
            
            % Partition system of equations
            Kfc = model.K(1:model.neqfree,model.neqfree+1:model.neq);
            Kcf = model.K(model.neqfree+1:model.neq,1:model.neqfree);
            Kcc = model.K(model.neqfree+1:model.neq,model.neqfree+1:model.neq);
            Ff  = model.F(1:model.neqfree);
            Dc  = model.D(model.neqfree+1:model.neq);
            Fc  = model.F(model.neqfree+1:model.neq);
            
            % Solve for Df
            Df = Kff \ (Ff - Kfc * Dc);
            
            % Reconstruct the global unknown vector D
            model.D = [ Df
                        Dc ];
            
            % Recover forcing unknown values (reactions) at essential B.C.
            % It is assumed that the Fc vector currently stores combined nodal
            % loads applyed directly to fixed d.o.f's.
            % Superimpose computed reaction values to combined nodal loads,
            % with inversed direction, that were applyed directly to fixed d.o.f.'s.
            Fc = -Fc + Kcf * Df + Kcc * Dc;
            
            % Reconstruct the global forcing vector
            model.F = [ Ff
                        Fc ];
        end
        
        %------------------------------------------------------------------
        % Processes current model data.
        function process(model)
            fprintf(1,'Preparing analysis data...\n');
            
            % Dimension and initialize global matrix and vectors
            model.dimKFD();
            
            % Generate global d.o.f. numbering matrix
            model.anm.setupDOFNum(model);
            
            % Assemble global d.o.f. numbering matrix
            model.assembleDOFNum();
            
            % Assemble element gather vector
            % (stores element d.o.f. eqn. numbers)
            model.assembleGle();
            
            % Store prescribed displacements in global displacement vector
            model.anm.setupPrescDispl(model);
            
            % Assemble global stiffness matrix
            fprintf(1,'Assembling stiffness matrix...\n');
            model.gblStiffMtx();
            
            % Assemble global forcing vector
            fprintf(1,'Assembling forcing vector...\n');
            model.anm.nodalLoads(model); % Initialize forcing vector with nodal loads
            model.edgeLoads();           % add edge (element side) loads to forcing vector
            model.areaLoads();           % add area (element) loads to forcing vector
            
            % Partition and solve system of equations
            fprintf(1,'Solving system of equations...\n');
            stability = model.solveEqnSystem();
            
            if (stability == 1)
                % Add element internal forces to global arrays of
                % internal forces
                fprintf(1,'Computing element internal forces...\n');
                %model.elemIntForce();

                % Add element displacements to array of internal displacements
                fprintf(1,'Computing element internal displacements...\n');
                %model.elemIntDispl();
            else
                fprintf(1,'Unstable structure or invalid input data...\n');
            end
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Model object.
        function clean(model)
            model.anm = [];
            model.elems = [];
            model.nodes = [];
            model.materials = [];
            model.nel = 0;
            model.nnp = 0;
            model.nmat = 0;
            model.neq = 0;
            model.neqfree = 0;
            model.neqfixed = 0;
            model.ID = [];
            model.K = [];
            model.F = [];
            model.D = [];
        end
    end
end