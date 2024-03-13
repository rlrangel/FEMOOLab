%% Anm_PlaneStrain class
classdef Anm_PlaneStrain < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_PlaneStrain()
            include_gblrefs;
            anm = anm@Anm(PLANE_STRAIN,2);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <anm.html *Anm*>.
    methods
        %------------------------------------------------------------------
        % Sets up global d.o.f (degree of freedom) numbering ID matrix.
        % Initializes ID matrix:
        %  if ID(k,n) = 0, d.o.f. k of node n is free.
        %  if ID(k,n) = 1, d.o.f. k of node n is constrained by support.
        % Counts total number of equations of free d.o.f.'s and total number
        % of equations of fixed d.o.f.'s.
        % Input arguments:
        %  model: handle to an object of the Model class
        function setupDOFNum(anm,model)
            % Dimension global d.o.f. numbering matrix
            model.ID = zeros(anm.ndof,model.nnp);
            
            % Compute total number of fixed d.o.f.
            model.neqfixed = 0;
            
            for n = 1:model.nnp
                if model.nodes(n).ebc(1) ~= 0
                    model.neqfixed = model.neqfixed + 1;
                    model.ID(1,n) = 1;
                end
                
                if model.nodes(n).ebc(2) ~= 0
                    model.neqfixed = model.neqfixed + 1;
                    model.ID(2,n) = 1;
                end
            end
            
            % Compute total number of free d.o.f.
            model.neqfree = model.neq - model.neqfixed;
        end
        
        %------------------------------------------------------------------
        % Stores prescribed displacements (known support settlement values)
        % in global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of the global displacement vector that corresponds to a free d.o.f.
        % Input arguments:
        %  model: handle to an object of the model class
        function setupPrescDispl(~,model)
            for n = 1:model.nnp
                if size(model.nodes(n).prescDispl,2) > 0
                    if (model.ID(1,n) > model.neqfree) && (model.nodes(n).prescDispl(1) ~= 0)
                        model.D(model.ID(1,n)) = model.nodes(n).prescDispl(1);
                    end
                    
                    if (model.ID(2,n) > model.neqfree) && (model.nodes(n).prescDispl(2) ~= 0)
                        model.D(model.ID(2,n)) = model.nodes(n).prescDispl(2);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Generates material constituive matrix for a given element.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        % Output:
        %  E: material constituive matrix
        function E = EMtx(~,elem)
            v = elem.material.poisson;
            elast = elem.material.elasticity;
            elast_coef = elast / ((1 + v) * (1 - (2*v)));
            
            E = elast_coef * [ 1-v  v    0;
                               v    1-v  0;
                               0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Generates cartesian nodal coordinates matrix for a given element.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        % Output:
        %  X: cartesian nodal coordinates matrix
        function X = elemNodalCoord(~,elem)
            X = zeros(elem.nen,2);
            
            for i = 1:elem.nen
                X(i,:) = [ elem.nodes(i).coords(1) elem.nodes(i).coords(2) ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates strain-displacement matrix at given position of an element.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  GradNcar: shape functions derivatives in cartesian coordinates
        % Output arguments:
        %  B: strain-displacement matrix evaluated at given position
        function B = BMtx(anm,elem,GradNcar,~,~,~)
            B = zeros(3,anm.ndof*elem.nen);
            
            for i = 1:elem.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0;
                B(2,2*i-1) = 0;               B(2,2*i) = GradNcar(2,i);
                B(3,2*i-1) = GradNcar(2,i);   B(3,2*i) = GradNcar(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Adds nodal load components to the global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  model: handle to an object of the Model class
        function nodalLoads(~,model)
            for n = 1:model.nnp
                id = model.ID(1,n);
                model.F(id) = model.F(id) + model.nodes(n).nodalLoad(1);

                id = model.ID(2,n);
                model.F(id) = model.F(id) + model.nodes(n).nodalLoad(2);
            end
        end
    end
end