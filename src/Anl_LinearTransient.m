%% Anl_LinearTransient Class
%
%% Description
%
% This is a sub-class of the <Anl.html Anl> class for the implementation
% of *Linear Transient Analysis*.
%
classdef Anl_LinearTransient < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        scheme    Scheme = Scheme.empty;  % handle to object of Scheme class
        max_step  int32  = int32.empty;   % maximum number of steps
        time_incr double = double.empty;  % time increment in each step
    end
    
    %% Constructor method
    methods
        function this = Anl_LinearTransient()
            this = this@Anl(Anl.LINEAR_TRANSIENT);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function status = Process(this,mdl)
            status = 1;
            
            % Check for any free d.o.f
            if (mdl.neqf == 0)
                fprintf('Model with no free degree-of-freedom!\n');
                status = 0;
                return;
            end
            
            % Assemble global matrices
            if (echo > 0)
                fprintf('Assembling global arrays...\n');
            end
            K = mdl.anm.GblStiffMtx(mdl); % stiffness matrix
            C = mdl.anm.GblRate1Mtx(mdl); % matrix related to 1st time derivatives of state variables
            M = mdl.anm.GblRate2Mtx(mdl); % matrix related to 2nd time derivatives of state variables
            
            % Assemble global forcing vector (currently assumed constant on time)
            F = zeros(mdl.neq,1);
            F = mdl.AddPointForce(F);
            F = mdl.AddEquivForce(F);
            
            % Assemble vector of fixed d.o.f.'s
            Uc = zeros(mdl.neq,1);
            Uc = mdl.AddPrescDOF(Uc);
            Uc = Uc(mdl.neqf+1:end);
            
            % Extract free d.o.f.'s of global arrays
            Kff = K(1:mdl.neqf,1:mdl.neqf);
            Cff = C(1:mdl.neqf,1:mdl.neqf);
            Mff = M(1:mdl.neqf,1:mdl.neqf);
            Ff  = F(1:mdl.neqf,1);
            
            % Apply EBCs by static condensation of forcing vector
            Kfc = K(1:mdl.neqf,mdl.neqf+1:end);
            Ff  = Ff - Kfc * Uc;
            
            % Assemble global initial conditions matrix
            IC = mdl.GblInitCondMtx();
            
            % Solve transient problem
            if (echo > 0)
                fprintf('Solving system of differential equations...\n');
            end
            dt   = this.time_incr;
            endt = this.time_incr * this.max_step;
            row  = mdl.neqf;
            col  = this.max_step+1;
            [U,Ut,Utt,steps,times] = this.scheme.Execute(dt,endt,row,col,IC,Kff,Cff,Mff,Ff);
            fprintf('Analysis completed!\n');
            
            % Initialize results vectors of state variables and time derivatives
            mdl.res.U   = zeros(mdl.neq,steps+1);
            mdl.res.Ut  = zeros(mdl.neq,steps+1);
            mdl.res.Utt = zeros(mdl.neq,steps+1);
            
            % Store results
            mdl.res.steps        = steps;
            mdl.res.times        = times;
            mdl.res.U(1:row,:)   = U;
            mdl.res.Ut(1:row,:)  = Ut;
            mdl.res.Utt(1:row,:) = Utt;
            
            % Add prescribed values of fixed d.o.f.'s
            mdl.res.U(row+1:end,:) = repmat(Uc,1,steps+1);
            
            % Zero out time derivatives of fixed d.o.f.'s
            mdl.res.Ut(row+1:end,:)  = zeros(mdl.neqc,steps+1);
            mdl.res.Utt(row+1:end,:) = zeros(mdl.neqc,steps+1);
        end
        
        %------------------------------------------------------------------
        function PosProcess(~,~)
            % Currently, no derived variables are calculated for
            % transient analysis, only state variables results are provided
            return;
        end
    end
end