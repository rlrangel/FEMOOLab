%% Simulation Class
%
%% Description
%
% This class defines a simulation object in the StAnOOP program.
% A simulation object is the one with the highest hierarchical level.
% Its properties are the three fundamental objects to perform a Finite
% Element analysis: <model.html Model>, <anl.html Analysis>, and
% <result.html Result>.
%
classdef Simulation < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        mdl = [];           % object of Model class
        anl = [];           % object of Anl (analysis) class
        res = drv.Result(); % object of Result (plot results driver) class
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function sim = Simulation()
            return;
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Open files and execute each simulation.
        function runAll(sim)
            % Get input file names
            [file,path] = uigetfile('*.*','StAnOOP - Input file','MultiSelect','on');
            if (isequal(file,0))
                return;
            end
            if (~iscell(file))
                file = {file};
            end
            
            % Open input files
            nfiles = length(file);
            fin = zeros(nfiles,1);
            for i=1:nfiles
                namein = strcat(path,string(file(i)));
                fin(i) = fopen(namein,'rt');
                if (fin(i) < 0)
                    fprintf('Error opening input file: %s\n',namein);
                    return;
                end
            end
            
            % Execute simulations
            for i=1:nfiles
                fprintf('Start of simulation %d: %s\n',i,string(file(i)));
                if (sim.runSim(fin(i)) == 0)
                    return;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Run one simulation: read file, pre-process/process/pos-process data,
        % and print/plot results.
        function status = runSim(sim,fin)
            status = 1;
            tic
            
            % Read input file
            fprintf('Reading model information...\n');
            read = drv.Read();
            if (read.inputFile(fin,sim) == 0)
                return;
            end
            
            % Pre-process
            fprintf('Preparing analysis data...\n');
            sim.anl.preProcess(sim.mdl);
            
            % Perform analysis
            if (sim.anl.process(sim.mdl,sim.res) == 0)
                return;
            end
            
            % Pos-process
            fprintf('Computing stresses...\n');
            sim.anl.posProcess(sim.mdl,sim.res);
            
            % Print results (TO DO)
%             fprintf('Analysis finished!Time: %.5f\n\n',toc);
%             nameout = strcat(namein(1:end-4),'.pos');
%             fout = fopen(nameout,'w');
%             %sim.res.print(mdl,fout);
%             fclose(fout);

            % Plot results
            fprintf('Plotting results...\n');
            sim.res.plot(sim.mdl);
            fprintf('Finished!\n');
        end
    end
end