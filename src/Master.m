%% Master class
%
%% Description
%
% This is the main class for running simulations and tests in FEMOOLab.
%
% It is responsible for managing the high-level tasks and call the
% appropriate methods to perform each stage of a simulation,
% from the reading of input files to the showing of results.
%
classdef Master < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        echo     int32   = int32.empty;    % echo level (amount of information to be printed in command window)
        n_files  int32   = int32.empty;    % number of input files
        path_in  string  = string.empty;   % path to input files folder
        files    string  = string.empty;   % full name of input files (with path)
        cur_file string  = string.empty;   % name of current file being run
        is_last  logical = logical.empty;  % flag for last input file to run
    end

    %% Constructor method
    methods
        function this = Master()
            this.SetDefaultProps();
        end
    end
    
    %% Public methods: main function
    methods
        %------------------------------------------------------------------
        function SetDefaultProps(this)
            this.echo = 1;
        end

        %------------------------------------------------------------------
        function RunSimulations(this,echo)
            if (nargin > 1)
                this.echo = echo;
            end

            % Print header
            this.PrintHeader();

            % Get input files
            if (~this.GetInputFiles())
                return;
            end

            % Display input files
            this.DisplayInputFiles();

            % Run each input file
            for i = 1:this.n_files
                this.is_last = (i == this.n_files);

                % Get current file name and extension
                this.cur_file = this.files(i);
                [~,~,ext] = fileparts(this.cur_file);

                % Run simulation
                if (strcmp(ext,'.json'))
                    [status,drv] = this.RunAnalysis();
                    if (~status)
                        continue;
                    end
                else
                    fprintf('Input file:\n%s\n',this.cur_file);
                    fprintf(2,'\nInvalid input file extension.\n');
                    this.PrintExit();
                    continue;
                end

                % Pos-process
                drv.PosProcess();
                
                % Finish current analysis
                fprintf('\nFinished!\n');
                if (~this.is_last)
                    fprintf('\n------------------------------------------------------------------\n\n');
                end
            end
        end

        %------------------------------------------------------------------
        function runTests(this,mode)
            this.echo = 0;
            
            % Print header
            this.PrintHeader();
            
            % Check input
            if (mode ~= 1 && mode ~= 2)
                fprintf('Invalid testing mode!\n');
                fprintf('\nExiting program...\n');
                return;
            end
            
            % Get input files
            if (~this.GetInputFiles())
                return;
            end
            
            % Display input files
            this.DisplayInputFiles();
            
            % Run each input file
            for i = 1:this.n_files
                this.is_last = (i == this.n_files);
                
                % Get current file name and extension
                this.cur_file = this.files(i);
                [~,~,ext] = fileparts(this.cur_file);
                if (~strcmp(ext,'.json'))
                    fprintf('Input file:\n%s\n',this.cur_file);
                    fprintf(2,'\nInvalid input file extension.\n');
                    this.PrintExit();
                    continue;
                end
                
                % Run simulation
                [status,drv] = this.RunAnalysis();
                if (~status)
                    continue;
                end
                
                % Pos-process
                drv.PosProcess();
                
                % Perform action according to testing mode
                if (mode == 1) % compare results file with reference
                    this.CompareTest(drv);
                elseif (mode == 2) % generate/update reference results
                   this.RenameRefResultFile(drv); 
                end
                
                % Finish current analysis
                if (~this.is_last)
                    fprintf('\n------------------------------------------------------------------\n\n');
                end
            end
            fprintf('\nFinished!\n');
        end

        %------------------------------------------------------------------
        function [status,drv] = RunAnalysis(this)
            % Open parameters file
            fprintf('Parameters file:\n%s\n',this.cur_file);
            fid = fopen(this.cur_file,'rt');
            if (fid < 0)
                fprintf(2,'\nError opening parameters file.\n');
                this.PrintExit();
                status = 0;
                return;
            end
            
            % Read parameters file
            if (this.echo > 0)
                fprintf('\nReading parameters file...\n');
            end
            read = Read();
            [status,drv] = read.Execute(this.path_in,fid);
            if (status == 0)
                this.PrintExit();
                return;
            end

            % Check input data
            if (this.echo > 0)
                fprintf('\nChecking consistency of input data...\n');
            end
            if (~read.check(drv))
                this.PrintExit();
                status = 0;
                return;
            end

            % Pre-process
            if (this.echo > 0)
                fprintf('\nPre-processing...\n');
            end
            if (~drv.preProcess())
                this.PrintExit();
                status = 0;
                return;
            end

            % Print simulation information
            if (this.echo > 0)
                this.PrintSimulationInfo(drv);
                fprintf('\nStarting analysis:\n');
                fprintf('%s\n',datestr(now));
            end
            
            % Show starting configuration
            if (this.echo > 0)
                Animation().CurConfig(drv,'Starting');
            end
            
            % Execute analysis
            tic;
            drv.Process();
            
            % Print finished status
            if (this.echo > 0)
                this.PrintFinishedStatus(drv,status);
            end
        end
    end
    
    %% Public methods: auxiliary functions
    methods
        %------------------------------------------------------------------
        function PrintHeader(~)
            fprintf('===================================================================\n');
            fprintf('    FEMOOLab - Finite Element Method Object-Oriented Laboratory    \n');
            fprintf('                    Version 1.0 - January 2023                     \n');
            fprintf('===================================================================\n\n');
        end
        
        %------------------------------------------------------------------
        function status = GetInputFiles(this)
            status = 1;
            
            % Get files from dialog
            filter  = {'*.json','Parameters File (*.json)'};
            title   = 'FEMOOLab - Input file';
            default = 'ProjectParameters.json';
            [file_names,path] = uigetfile(filter,title,default,'MultiSelect','on');
            if (isequal(file_names,0))
                fprintf('No file selected.\n');
                fprintf('\nExiting program...\n');
                status = 0;
                return;
            end
            file_fullnames = fullfile(path,file_names);
            
            % Convert to string array
            this.files   = string(file_fullnames);
            this.path_in = string(path);
        end

        %------------------------------------------------------------------
        function DisplayInputFiles(this)
            this.n_files = length(this.files);
            fprintf('%d input files selected:\n',this.n_files);
            for i = 1:this.n_files
                fprintf('%s\n',this.files(i));
            end
            fprintf('\n------------------------------------------------------------------\n\n');
        end

        %------------------------------------------------------------------
        function PrintSimulationInfo(~,drv)
            fprintf('\nSimulation ready:\n');
            if (~isempty(drv.name))
                fprintf('Name...................: %s\n',drv.name);
            end

            % TODO
        end
        
        %------------------------------------------------------------------
        function PrintFinishedStatus(~,drv,status)
            fprintf('\n\nAnalysis finished:\n');
            fprintf('%s\n',datestr(now));
            if (status == 1)
                time = seconds(toc);
                time.Format = 'hh:mm:ss.SS';
                fprintf('Total time: %s\n',string(time));
            elseif (status == 2)
                curr_time  = seconds(toc);
                total_time = seconds(drv.total_time);
                curr_time.Format  = 'hh:mm:ss.SS';
                total_time.Format = 'hh:mm:ss.SS';
                fprintf('Current analysis time: %s\n',string(curr_time));
                fprintf('Total simulation time: %s\n',string(total_time));
            end
        end
        
        %------------------------------------------------------------------
        function PrintExit(this)
            if (this.is_last)
                fprintf('\nExiting program...\n');
            else
                fprintf('\nAborted!\n');
                fprintf('\n------------------------------------------------------------------\n\n');
            end
        end

        %------------------------------------------------------------------
        function CompareTest(this,drv)
            % Check if reference and current results files exist
            name_ref = strcat(drv.path_in,drv.name,"_ref.pos");
            name_cur = strcat(drv.path_out,drv.name,".pos");
            if (exist(name_ref,'file') ~= 2)
                fprintf(2,'\nMissing reference results file!\n');
                this.DeleteFile(name_cur);
                status = rmdir(drv.path_out); %#ok<NASGU>
                return;
            end
            if (exist(name_cur,'file') ~= 2)
                fprintf(2,'\nCurrent results file was not generated correctly!\n');
                this.DeleteFile(name_cur);
                status = rmdir(drv.path_out); %#ok<NASGU>
                return;
            end
            
            % Compare contents of files
            file_ref = javaObject('java.io.File',name_ref);
            file_cur = javaObject('java.io.File',name_cur);
            is_equal = javaMethod('contentEquals','org.apache.commons.io.FileUtils',file_ref,file_cur);
            if (is_equal)
                fprintf(1,'\nTest passed!\n');
            else
                fprintf(2,'\nResults are different!\n');
            end
            
            % Delete current results file and folder (if empty)
            this.DeleteFile(name_cur);
            status = rmdir(drv.path_out); %#ok<NASGU>
        end
        
        %------------------------------------------------------------------
        function RenameRefResultFile(this,drv)
            % Check if current results file exist
            name_cur = strcat(drv.path_out,drv.name,".pos");
            if (exist(name_cur,'file') ~= 2)
                fprintf(2,'\nCurrent results file was not generated correctly!\n');
                this.DeleteFile(name_cur);
                status = rmdir(drv.path_out); %#ok<NASGU>
                return;
            end
            
            % Move current results file out of output folder
            if (~movefile(name_cur,drv.path_in))
                fprintf(2,'\nCurrent results file was not generated correctly!\n');
                this.DeleteFile(name_cur);
                status = rmdir(drv.path_out); %#ok<NASGU>
                return;
            end
            
            % Delete output folder (if empty)
            status = rmdir(drv.path_out); %#ok<NASGU>
            
            % Rename current results file to reference results file
            name_cur = strcat(drv.path_in,drv.name,".pos");
            name_ref = strcat(drv.path_in,drv.name,"_ref.pos");
            movefile(name_cur,name_ref);
            fprintf(1,'\nReference results generated!\n');
        end
        
        %------------------------------------------------------------------
        function DeleteFile(~,file_name)
            if (exist(file_name,'file') == 2)
                warning off MATLAB:DELETE:FileNotFound
                delete(sprintf('%s',file_name));
                warning on MATLAB:DELETE:FileNotFound
            end
        end
    end
end