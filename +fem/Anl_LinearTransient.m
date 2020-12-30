%% Anl_LinearTransient Class (Linear Transient Analysis)
%
%% Description
%
% This is a sub-class in the FEMOOLab program that implements abstract 
% methods declared in <anl.html Anl: analysis super-class> to deal
% with linear-transient analysis.
%
%% Class definition
%
classdef Anl_LinearTransient < fem.Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        scheme = [];                     % object of Scheme class
        incr     double = double.empty;  % increment of time in each step
        max_step int32  = int32.empty;   % maximum number of steps
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_LinearTransient()
            this = this@fem.Anl(fem.Anl.LINEAR_TRANSIENT);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process linear-transient analysis by assembling global semi-
        % discretized system of ODEs on time and calculating state variables
        % and their time derivatives.
        function status = process(this,sim)
            status = 1;
            mdl = sim.mdl;
            res = mdl.res;
            
            % Check for any free d.o.f
            if (mdl.neqf == 0)
                fprintf('Model with no free degree-of-freedom!\n');
                status = 0;
                return;
            end
            
            % Assemble global matrices
            fprintf('Assembling global arrays...\n');
            K = mdl.anm.gblStiffMtx(mdl); % stiffness matrix
            C = mdl.anm.gblVelMtx(mdl);   % matrix related to 1st time derivative ("velocity" matrix)
            M = mdl.anm.gblAccelMtx(mdl); % matrix related to 2nd time derivative ("acceleration" matrix)
            
            % Assemble global forcing vector (currently assumed constant on time)
            F = zeros(mdl.neq,1);
            F = mdl.anm.addPointForce(mdl,F);
            F = mdl.anm.addEquivForce(mdl,F);
            
            % Extract free-free terms of global arrays
            K = K(1:mdl.neqf,1:mdl.neqf);
            C = C(1:mdl.neqf,1:mdl.neqf);
            M = M(1:mdl.neqf,1:mdl.neqf);
            F = F(1:mdl.neqf,1);
            
            % Assemble global initial conditions matrix
            IC = mdl.anm.gblInitCondMtx(mdl);
            
            % Prepare input of scheme method
            dt   = this.incr;
            endt = this.incr * this.max_step;
            row  = mdl.neqf;
            col  = this.max_step+1;
            
            % Solve transient problem
            fprintf('Solving system of differential equations...\n');
            [U,Ut,Utt,steps,times] = this.scheme.execute(dt,endt,row,col,IC,K,C,M,F);
            fprintf('Analysis completed!\n');
            
            % Initialize results vectors of state variables and time derivatives
            res.U   = zeros(mdl.neq,steps+1);
            res.Ut  = zeros(mdl.neq,steps+1);
            res.Utt = zeros(mdl.neq,steps+1);
            
            % Store results
            res.steps        = steps;
            res.times        = times;
            res.U(1:row,:)   = U;
            res.Ut(1:row,:)  = Ut;
            res.Utt(1:row,:) = Utt;
            
            % Add prescribed values of fixed d.o.f.'s
            Uc = zeros(mdl.neq,1);
            mdl.anm.addEBC(mdl,Uc);
            res.U(row+1:end,:) = repmat(Uc(row+1:end),1,steps+1);
            
            % Zero out time derivatives of fixed d.o.f.'s
            res.Ut(row+1:end,:)  = zeros(mdl.neqc,steps+1);
            res.Utt(row+1:end,:) = zeros(mdl.neqc,steps+1);
        end
    end
end