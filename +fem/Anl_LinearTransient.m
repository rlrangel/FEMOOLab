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
        max_step int32  = int32.empty;   % maximum number of steps
        incr     double = double.empty;  % increment of time in each step
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
        % Process linear transient analysis by assembling the global semi-
        % discretized system of ODEs on time and calculating state variables
        % and their time derivatives.
        function status = process(this,mdl)
            status = 1;
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
            C = mdl.anm.gblRate1Mtx(mdl); % matrix related to 1st time derivatives
            M = mdl.anm.gblRate2Mtx(mdl); % matrix related to 2nd time derivatives
            
            % Assemble global forcing vector (currently assumed constant on time)
            F = zeros(mdl.neq,1);
            F = mdl.addPointForce(F);
            F = mdl.addEquivForce(F);
            
            % Assemble vector of fixed d.o.f.'s
            Uc = zeros(mdl.neq,1);
            Uc = mdl.addPrescDOF(Uc);
            Uc = Uc(mdl.neqf+1:end);
            
            % Extract free d.o.f.'s terms of global arrays
            Kff = K(1:mdl.neqf,1:mdl.neqf);
            Cff = C(1:mdl.neqf,1:mdl.neqf);
            Mff = M(1:mdl.neqf,1:mdl.neqf);
            Ff  = F(1:mdl.neqf,1);
            
            % Static condensation of forcing vector (apply EBCs)
            Kfc = K(1:mdl.neqf,mdl.neqf+1:end);
            Ff  = Ff - Kfc * Uc;
            
            % Assemble global initial conditions matrix
            IC = mdl.gblInitCondMtx();
            
            % Solve transient problem
            fprintf('Solving system of differential equations...\n');
            dt   = this.incr;
            endt = this.incr * this.max_step;
            row  = mdl.neqf;
            col  = this.max_step+1;
            [U,Ut,Utt,steps,times] = this.scheme.execute(dt,endt,row,col,IC,Kff,Cff,Mff,Ff);
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
            res.U(row+1:end,:) = repmat(Uc,1,steps+1);
            
            % Zero out time derivatives of fixed d.o.f.'s
            res.Ut(row+1:end,:)  = zeros(mdl.neqc,steps+1);
            res.Utt(row+1:end,:) = zeros(mdl.neqc,steps+1);
        end
        
        %------------------------------------------------------------------
        % Pos-process results to compute derived variables.
        function posProcess(~,~)
            % Currently, no derived variables are calculated for
            % transient analysis, only state variables results are provided
            return;
        end
    end
end