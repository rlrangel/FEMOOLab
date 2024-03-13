% Clear memory
clc
clear
close all

% Create model object
model = Model();

% Open input neutral-format file
filterspec = {'*.nf';'*.dat';'*.pos'};
[filename,pathname] = uigetfile(filterspec,'Elasticity2D - Input file');

% Check for valid input file
if filename ~= 0
    fullname = strcat(pathname,filename);
    fid = fopen(fullname,'rt');
    
    % Read input FE model data
    fprintf(1,'Pre-processing...\n');
    preReadNF(fid,model);
    
    % Process provided FE model data
    model.process();
end
