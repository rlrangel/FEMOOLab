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
        function status = process(this,sim)
            status = 1;
            mdl = sim.mdl;
            res = mdl.res;
            
            % Assemble global stiffness matrix
            fprintf('Assembling stiffness matrix...\n');
            K = mdl.anm.gblStiffMtx(mdl);
            
            % Check model stability
            if (this.singularMtx(mdl,K))
                fprintf(1,'Singular stiffness matrix!\n');
                status = 0;
                return;
            end
            
            % Assemble global forcing vector
            fprintf('Assembling forcing vector...\n');
            F = zeros(mdl.neq,1);
            F = mdl.anm.addPointForce(mdl,F);
            F = mdl.anm.addEquivForce(mdl,F);
            
            % Assemble global vector of state variables
            U = zeros(mdl.neq,1);
            U = mdl.anm.addEBC(mdl,U);
            
            % Solve system of equations and store result
            fprintf('Solving system of equations...\n');
            res.U = this.solveSystem(mdl,K,F,U);
        end
    end
end