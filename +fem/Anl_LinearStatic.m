%% Anl_LinearStatic Class (Linear Estatic Analysis)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Anl: analysis super-class> to deal
% with linear-static analysis.
%
%% Class definition
%
classdef Anl_LinearStatic < fem.Anl
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_LinearStatic()
            this = this@fem.Anl(fem.Anl.LINEAR_STATIC);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process linear-static analysis by assembling global equilibrium
        % system of equations and calculating state variables.
        function status = process(this,mdl)
            status = 1;
            res = mdl.res;
            
            % Assemble global stiffness matrix
            fprintf('Assembling stiffness matrix...\n');
            K = mdl.gblStiffMtx();
            
            % Check for singular matrix
            if ((rcond(K(1:mdl.neqf,1:mdl.neqf)) < 10e-15))
                fprintf(1,'Singular stiffness matrix!\n');
                status = 0;
                return;
            end
            
            % Assemble global forcing vector
            fprintf('Assembling forcing vector...\n');
            F = zeros(mdl.neq,1);
            F = mdl.addPointForce(F);
            F = mdl.addEquivForce(F);
            
            % Assemble global vector of state variables
            U = zeros(mdl.neq,1);
            U = mdl.addPrescDOF(U);
            
            % Solve system of equations and store result
            fprintf('Solving system of equations...\n');
            res.U = this.solveSystem(mdl,K,F,U);
        end
        
        %------------------------------------------------------------------
        % Pos-process results to compute derived variables.
        function posProcess(~,mdl)
            % Set element transformation matrices of gauss-to-node results
            mdl.setupTGNmtx();
            
            % Initialize result arrays
            mdl.res.initPosResults(mdl);
            
            % Compute derived variables and principal values and directions
            % at Gauss points
            mdl.anm.gaussDerivedVar(mdl);
            
            % Extrapolate Gauss point results to element node results
            mdl.anm.elemDerivedVarExtrap(mdl);
            
            % Smooth element node result to global node results
            mdl.anm.nodeDerivedVarExtrap(mdl);
            
            % Compute minimum and maximum values of obtained results
            mdl.res.setMinMaxValues();
            
            % Clear numerical garbage
            mdl.res.clearSmallValues(mdl);
        end
    end
end