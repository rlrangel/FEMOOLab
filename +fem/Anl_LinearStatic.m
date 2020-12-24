%% Anl_LinearStatic Class
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Anl: analysis super-class> to deal
% with linear-static analysis.
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
        % Process linear-static analysis based on the direct stiffness
        % method, i.e., assemble global equilibrium system of equations,
        % and calculate state variables (nodal displacements).
        function status = process(this,mdl,res)
            status = 1;
            
            % Assemble global elastic stiffness matrix
            fprintf('Assembling stiffness matrix...\n');
            K = mdl.gblElastStiffMtx();
            
            % Check model stability
            if (this.singularMtx(mdl,K))
                fprintf(1,'Unstable model!\n');
                status = 0;
                return;
            end
            
            % Assemble global forcing vector
            fprintf('Assembling forcing vector...\n');
            F = zeros(mdl.neq,1);
            F = mdl.anm.addPointLoad(mdl,F);
            F = mdl.addEquivLoad(F);
            
            % Add prescribed displacements to global displacement vector
            D = zeros(mdl.neq,1);
            D = mdl.anm.addPrescDispl(mdl,D);
            
            % Solve system of equations and store result
            fprintf('Solving system of equations...\n');
            res.D = this.solveSystem(mdl,K,F,D);
        end
    end
end