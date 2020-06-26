%% Anm_Axisymmetric Class
%
% This is a sub-class in the StAnOOP program that implements abstract 
% methods declared in super-class Anm to deal with axisymmetric models.
%
classdef Anm_Axisymmetric < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Axisymmetric()
            c = Constants();
            anm = anm@Anm(c.AXISYMMETRIC,2);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anm
    methods
        %------------------------------------------------------------------
        % Assemble material constitutive matrix for a given element.
        function C = Cmtx(~,elem)
            E = elem.mat.E;
            v = elem.mat.v;
            e = E/((1+v)*(1-(2*v)));
            
            C = e * [ 1-v  v    v    0;
                      v    1-v  v    0;
                      v    v    1-v  0;
                      0    0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble strain-displacement matrix at a given position of
        % an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        %  r,s: parametric coordinates
        function B = Bmtx(anm,elem,GradNcar,r,s)
            % Evaluate matrix of shape functions
            N = elem.Nmtx(r,s);
            
            % Location of evaluation point (X coordinate is the radius in axisymm.)
            p = N * elem.carCoord;
            radius = p(1);
            
            % Assemble strain-displacement matrix
            B = zeros(4,elem.nen*anm.ndof);
            
            for i = 1:elem.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0.0;
                B(2,2*i-1) = N(i)/radius;     B(2,2*i) = 0.0;
                B(3,2*i-1) = 0.0;             B(3,2*i) = GradNcar(2,i);
                B(4,2*i-1) = GradNcar(2,i);   B(4,2*i) = GradNcar(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        % Input:
        %  C:   constituive matrix
        %  B:   strain-displacement matrix
        %  d:   generalized displacements for all d.o.f.'s of element
        function str = pointStress(~,C,B,d)
            % Compute point stress components
            str_raw = C * B * d;
            
            % Skip tangential stress component
            str(1) = str_raw(1);
            str(2) = str_raw(3);
            str(3) = str_raw(4);
        end
    end
end