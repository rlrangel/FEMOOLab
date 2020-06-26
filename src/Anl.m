%% Anl (Analysis) Class
%
% This is an abstract super-class that generically specifies an analysis 
% type in the StAnOOP program.
%
% Essentially, this super-class declares abstract methods that are the
% functions that should be implemented in a derived sub-class that deals
% with a specific type of analysis.
%
classdef Anl < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = 0;   % flag for type of analysis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl(type)
            anl.type = type;
        end
    end
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Process model data to compute results.
        status = process(anl,mdl,res);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Pre-process model data to setup mesh information.
        function preProcess(~,mdl)
            % Initialize global d.o.f. numbering matrix
            mdl.setupDOFNum();
            
            % Assemble global d.o.f. numbering matrix
            mdl.assembleDOFNum();
            
            % Assemble element gather vectors
            mdl.assembleGle();
        end
        
        %------------------------------------------------------------------
        % Pos-process results to compute derived quantitites.
        function posProcess(anl,mdl,res)
            % Compute stresses at Gauss points and principal stresses
            mdl.gaussStress(anl,res);
            
            % Extrapolate Gauss point results to element node results
            mdl.elemStressExtrap(res);
            
            % Smooth element node result to global node results
            mdl.nodeStressExtrap(res);
            
            % Clear numerical garbage
            res.clearSmallValues(mdl);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Check matrix singularity by checking its reciprocal condition number.
        % A very low reciprocal condition number indicates that the matrix
        % is badly conditioned and may be singular.
        function singular = singularMtx(mdl,K)
            if (rcond(K(1:mdl.neqf,1:mdl.neqf)) < 10e-15)
                singular = 1;
            else
                singular = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Partition and solve a linear system of equations.
        %  f --> free d.o.f. (natural B.C. - unknown) 
        %  c --> constrained d.o.f. (essential B.C. - known) 
        %
        % [ Kff Kfc ] * [ Df ] = [ Ff ]
        % [ Kcf Kcc ]   [ Dc ] = [ Fc ]
        %
        function D = solveSystem(mdl,K,F,D)
            % Free and constrained d.o.f. terms
            f = 1:mdl.neqf;
            c = mdl.neqf+1:mdl.neq;
            
            % Partition system of equations
            Kff = K(f,f);
            Kfc = K(f,c);
            Ff  = F(f);
            Dc  = D(c);
            
            % Solve for Df
            Df = Kff \ (Ff - Kfc * Dc);
            
            % Reconstruct unknown vector D
            D = [ Df; Dc ];
        end
        
        %------------------------------------------------------------------
        % Compute principal stress components and orientation for a given
        % stress tensor.
        % Input:
        %  str:  stress tensor (sx, sy, txy) stored in a column vector.
        % Output:
        %  prc:    principal stress components (s1, s2, taumax) stored in
        %          a column vector.
        %  thetap: angle of normal of principal stress plane w.r.t x axis
        %          (angle is returned in radians from 0 to 180 degrees).
        function [prc,thetap] = princStress(str)
            sx  = str(1);
            sy  = str(2);
            txy = str(3);
            
            center = (sx+sy)/2;
            deltas = (sx-sy)/2;
            radius = sqrt((deltas^2) + (txy*txy));
            
            prc = zeros(3,1);
            prc(1) = center + radius;   % s1
            prc(2) = center - radius;   % s2
            prc(3) = radius;            % taumax
            
            if (abs(deltas) > 0.0)
                thetap = 0.5 * atan2(txy,deltas);
            elseif (txy > 0.0)
                thetap = pi / 4.0;
            elseif (txy < 0.0)
                thetap = -pi / 4.0;
            else
                thetap = 0.0;
            end
            
            if( thetap < 0.0 )
                thetap = pi + thetap;
            end
        end
    end
end