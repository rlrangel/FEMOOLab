%% Anl_LinearSteadyState Class
%
%% Description
%
% This is a sub-class of the <Anl.html Anl> class for the implementation
% of *Linear Steady State Analysis*.
%
classdef Anl_LinearSteadyState < Anl
    %% Constructor method
    methods
        function this = Anl_LinearSteadyState()
            this = this@Anl(Anl.LINEAR_STEADYSTATE);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function status = Process(this,mdl,echo)
            status = 1;

            % Assemble global stiffness matrix
            if (echo > 0)
                fprintf('Assembling stiffness matrix...\n');
            end
            K = mdl.anm.GblStiffMtx(mdl);
            
            % Check for singular matrix
            singular = rcond(K(1:mdl.neqf,1:mdl.neqf)) < 10e-15;
            if (singular)
                fprintf(1,'Singular stiffness matrix!\n');
                status = 0;
                return;
            end
            
            % Assemble global forcing vector
            if (echo > 0)
                fprintf('Assembling forcing vector...\n');
            end
            F = zeros(mdl.neq,1);
            F = mdl.AddPointForce(F);
            F = mdl.AddEquivForce(F);
            
            % Add stabilization of convective term
            [K,F] = mdl.anm.StabConvec(mdl,K,F);
            
            % Assemble global vector of state variables
            U = zeros(mdl.neq,1);
            U = mdl.AddPrescDOF(U);
            
            % Solve system of equations and store results
            if (echo > 0)
                fprintf('Solving system of equations...\n');
            end
            mdl.res.U = this.SolveSystem(mdl,K,F,U);
        end
        
        %------------------------------------------------------------------
        function PosProcess(~,mdl)
            % Initialize results arrays
            mdl.res.InitPosResults(mdl);
            
            % Compute derived variables and principal values at Gauss points
            mdl.anm.GaussDerivVar(mdl);
            
            % Extrapolate Gauss point results to element node results
            mdl.anm.ElemDerivVarExtrap(mdl);
            
            % Smooth element node result to global node results
            mdl.anm.NodeDerivVarExtrap(mdl);
            
            % Compute minimum and maximum values of results
            mdl.res.SetMinMaxValues(mdl);
            
            % Clear numerical garbage
            mdl.res.ClearSmallValues(mdl);
        end
    end
end