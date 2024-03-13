%% Anm class (Analysis Model Class)
classdef Anm < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        analysis_type = 0;   % type of analysis model
        ndof = 0;            % number of degrees-of-freedom per node
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm(type,n)
            anm.analysis_type = type;
            anm.ndof = n;
        end
    end
    
    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Sets up global d.o.f (degree of freedom) numbering ID matrix.
        % Initializes ID matrix:
        %  if ID(k,n) = 0, d.o.f. k of node n is free.
        %  if ID(k,n) = 1, d.o.f. k of node n is constrained by support.
        % Counts total number of equations of free d.o.f.'s and total number
        % of equations of fixed d.o.f.'s.
        % Input arguments:
        %  model: handle to an object of the Model class
        setupDOFNum(anm,model)
        
        %------------------------------------------------------------------
        % Stores prescribed displacements (known support settlement values)
        % in global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of the global displacement vector that corresponds to a free d.o.f.
        % Input arguments:
        %  model: handle to an object of the model class
        setupPrescDispl(~,model)
        
        %------------------------------------------------------------------
        % Generates material constituive matrix for a given element.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        % Output:
        %  E: material constituive matrix
        E = EMtx(~,elem)
        
        %------------------------------------------------------------------
        % Generates cartesian nodal coordinates matrix for a given element.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        % Output:
        %  X: cartesian nodal coordinates matrix
        X = elemNodalCoord(~,elem)
        
        %------------------------------------------------------------------
        % Generates strain-displacement matrix at given position of an element.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  GradNcar: shape functions derivatives in cartesian coordinates
        %  X: cartesian nodal coordinates matrix
        %  r,s: parametric coordinate values
        % Output arguments:
        %  B: strain-displacement matrix evaluated at given position
        B = BMtx(anm,elem,GradNcar,X,r,s)
        
        %------------------------------------------------------------------
        % Adds nodal load components to the global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  model: handle to an object of the Model class
        nodalLoads(~,model)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of an Anm object.
        function clean(anm)
            anm.analysis_type = 0;
            anm.ndof = 0;
        end
    end
end