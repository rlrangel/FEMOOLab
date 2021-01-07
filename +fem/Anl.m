%% Anl Class (Analysis)
%
%% Description
%
% This is an abstract super-class that generically specifies an analysis 
% type in the FEMOOLab program.
% Essentially, this super-class declares abstract methods that define
% the general behavior of analysis types. These abstract methods
% are the functions that should be implemented in a derived sub-class
% that deals with specific types of analysis.
%
%% Subclasses
%
% * <anl_linearstatic.html Anl_LinearStatic: linear static analysis subclass>
% * <anl_lineartransient.html Anl_LinearTransient: linear transient analysis subclass>
%
%% Class definition
%
classdef Anl < handle
    %% Constant values
    properties (Constant = true, Access = public)
        % Types of analysis
        LINEAR_STATIC    = int32(1);
        LINEAR_TRANSIENT = int32(2);
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
        % Process model data to compute state variables results.
        status = process(this,mdl);
        
        %------------------------------------------------------------------
        % Pos-process results to compute derived variables.
        posProcess(this,mdl);
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Pre-process model data to setup mesh information.
        function preProcess(~,mdl)
            % Compute total number of equations
            mdl.neq = mdl.nnp * mdl.anm.ndof;
            
            % Initialize global d.o.f. numbering matrix
            mdl.setupDOFNum();
            
            % Assemble global d.o.f. numbering matrix
            mdl.assembleDOFNum();
            
            % Assemble element gather vectors
            mdl.assembleGle();
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
    end
end