%% Anl_LinearElastic Class
%
%% Description
%
% This is a sub-class in the StAnOOP program that implements abstract 
% methods declared in <anl.html Anl: analysis super-class> to deal
% with linear-elastic analysis.
%
classdef Anl_LinearElastic < fem.Anl
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_LinearElastic()
            this = this@fem.Anl(fem.Anl.LINEAR_ELASTIC);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process linear-elastic analysis based on the direct stiffness
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