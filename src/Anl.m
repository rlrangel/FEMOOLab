%% Anl Class
%
%% Description
%
% This is a handle super-class for the definition of finite element
% analysis types.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <Anl_LinearSteadyState.html Anl_LinearSteadyState>
% * <Anl_LinearTransient.html Anl_LinearTransient>
%
classdef Anl < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of analysis
        LINEAR_STEADYSTATE = int32(1);
        LINEAR_TRANSIENT   = int32(2);
    end
    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type int32 = int32.empty;  % flag for type of analysis
    end
    
    %% Constructor method
    methods
        function this = Anl(type)
            if (nargin > 0)
                this.type = type;
            end
        end
    end
    
    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Assemble global system of equilibrium equations and calculate
        % state variables.
        status = Process(this,mdl);
        
        %------------------------------------------------------------------
        % Compute derived variables.
        PosProcess(this,mdl);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function PreProcess(~,mdl)
            % Compute total number of equations
            mdl.neq = mdl.nnp * mdl.anm.ndof;
            
            % Initialize global d.o.f. numbering matrix
            mdl.SetupDOFNum();
            
            % Assemble global d.o.f. numbering matrix
            mdl.AssembleDOFNum();
            
            % Assemble element gather vectors
            mdl.AssembleGle();
            
            % Set element initial properties
            mdl.anm.SetupElemProps(mdl);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Partition and solve a linear system of equations.
        %  f --> free d.o.f. (natural B.C. - unknown) 
        %  c --> constrained d.o.f. (essential B.C. - known) 
        %
        % [ Kff Kfc ] * [ Uf ] = [ Ff ]
        % [ Kcf Kcc ]   [ Uc ] = [ Fc ]
        %
        function U = SolveSystem(mdl,K,F,U)
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
    end
end