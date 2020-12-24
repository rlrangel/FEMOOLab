%% Simulation Class
%
%% Description
%
% This class defines a simulation object in the FEMOOLab program.
% A simulation object is the one with the highest hierarchical level.
% Its properties are the two fundamental objects to perform a Finite
% Element analysis: <model.html Model> and <anl.html Analysis>.
%
%% Class definition
%
classdef Simulation < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        mdl  = [];  % object of Model class
        anl  = [];  % object of Anl (Analysis) class
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Simulation(opt)
            this.runFiles(opt);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Open files and execute each simulation.
        function runFiles(sim,opt)
            % Get input file names
            [file_names,file_path] = uigetfile('*.*','FEMOOLab - Input file','MultiSelect','on');
            if (isequal(file_names,0))
                return;
            end
            if (~iscell(file_names))
                file_names = {file_names};
            end
            
            % Open and read input files to run each simulation
            for i = 1:length(file_names)
                % Get file parts name
                full_name = strcat(file_path,string(file_names(i)));
                [~,name,~] = fileparts(full_name);
                
                fprintf('Start of analysis: %s\n',name);
                
                % Get file id
                fid = fopen(full_name,'rt');
                if (fid < 0)
                    fprintf('Error opening input file: %s\n\n\n',name);
                    continue;
                end
                
                % Create Model object
                sim.mdl = fem.Model();
                
                % Read input file
                fprintf('Reading model information...\n');
                if (~drv.Read().execute(fid,sim,opt))
                    fprintf('\n\n');
                    continue;
                end
                
                % Execute simulation
                sim.run();
                fprintf('\n\n');
            end
        end
        
        %------------------------------------------------------------------
        % Run a simulation and plot/print results, assuming that the
        % simulation object has already been filled correctly.
        function run(sim)
            % Pre-process
            fprintf('Preparing analysis data...\n');
            sim.anl.preProcess(sim.mdl);
            
            % Perform analysis
            if (sim.anl.process(sim))
                % Pos-process
                fprintf('Computing stresses...\n');
                sim.anl.posProcess(sim.mdl);
                
                % Plot results
                fprintf('Plotting results...\n');
                sim.mdl.res.plot(sim.mdl);
                fprintf('Finished!\n');
                
                % Print results... (TO DO)
            end
        end
    end
end