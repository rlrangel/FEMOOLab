%% Anm_PlaneStrain Class
%
% This is a sub-class in the StAnOOP program that implements abstract 
% methods declared in super-class Anm to deal with plane strain models.
%
classdef Anm_PlaneStrain < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_PlaneStrain()
            c = Constants();
            anm = anm@Anm(c.PLANE_STRAIN,2);
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
            
            C = e * [ 1-v  v    0;
                      v    1-v  0;
                      0    0    (1-(2*v))/2 ];
        end
        
        %------------------------------------------------------------------
        % Assemble strain-displacement matrix at a given position of
        % an element.
        % Input:
        %  GradNcar: shape functions derivatives w.r.t. cartesian coordinates
        function B = Bmtx(anm,elem,GradNcar,~,~)
            B = zeros(3,elem.nen*anm.ndof);
            
            for i = 1:elem.nen
                B(1,2*i-1) = GradNcar(1,i);   B(1,2*i) = 0;
                B(2,2*i-1) = 0;               B(2,2*i) = GradNcar(2,i);
                B(3,2*i-1) = GradNcar(2,i);   B(3,2*i) = GradNcar(1,i);
            end
        end
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        % Input:
        %  C:   constituive matrix
        %  B:   strain-displacement matrix
        %  d:   generalized displacements for all d.o.f.'s of element
        function str = pointStress(~,C,B,d)
            % In plane strain, raw stress vector is the target one
            str = C * B * d;
        end
    end
end