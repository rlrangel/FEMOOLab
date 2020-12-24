%% Anl Class (Analysis)
%
%% Description
%
% This is an abstract super-class that generically specifies an analysis 
% type in the FEMOOLab program.
%
% Essentially, this super-class declares abstract methods that are the
% functions that should be implemented in a derived sub-class that deals
% with a specific type of analysis. These abstract methods
% are the functions that should be implemented in a derived sub-class
% that deals with specific types of analysis.
%
%% Subclasses
%
% * <anl_linearstatic.html Anl_LinearStatic: linear static analysis subclass>
%
%% Class definition
%
classdef Anl < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of analysis
        GENERIC       = int32(0);
        LINEAR_STATIC = int32(1);    
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;  % flag for type of analysis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl(type)
            this.type = type;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Process model data to compute results.
        status = process(anl,sim);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Pre-process model data to setup mesh information.
        function preProcess(~,mdl)
            % Compute total number of equations
            mdl.neq = mdl.nnp * mdl.anm.ndof;
            
            % Initialize global d.o.f. numbering matrix
            mdl.anm.setupDOFNum(mdl);
            
            % Assemble global d.o.f. numbering matrix
            mdl.assembleDOFNum();
            
            % Assemble element gather vectors
            mdl.assembleGle();
        end
        
        %------------------------------------------------------------------
        % Pos-process results to compute derived quantitites.
        function posProcess(anl,mdl)
            % Compute stresses at Gauss points and principal stresses
            mdl.gaussStress(anl);
            
            % Extrapolate Gauss point results to element node results
            mdl.elemStressExtrap();
            
            % Smooth element node result to global node results
            mdl.nodeStressExtrap();
            
            % Clear numerical garbage
            mdl.res.clearSmallValues(mdl);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Check matrix singularity by checking its reciprocal condition number.
        % A very low reciprocal condition number indicates that the matrix
        % is badly conditioned and may be singular.
        function singular = singularMtx(mdl,K)
            singular = (rcond(K(1:mdl.neqf,1:mdl.neqf)) < 10e-15);
        end
        
        %------------------------------------------------------------------
        % Partition and solve a linear system of equations.
        %  f --> free d.o.f. (natural B.C. - unknown) 
        %  c --> constrained d.o.f. (essential B.C. - known) 
        %
        % [ Kff Kfc ] * [ Uf ] = [ Ff ]
        % [ Kcf Kcc ]   [ Uc ] = [ Fc ]
        %
        function U = solveSystem(mdl,K,F,U)
            % Free and constrained d.o.f. terms
            f = 1:mdl.neqf;
            c = mdl.neqf+1:mdl.neq;
            
            % Partition system of equations
            Kff = K(f,f);
            Kfc = K(f,c);
            Ff  = F(f);
            Uc  = U(c);
            
            % Solve for Uf
            Uf = Kff \ (Ff - Kfc * Uc);
            
            % Reconstruct unknown vector U
            U = [ Uf; Uc ];
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