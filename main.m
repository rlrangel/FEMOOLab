%% StAnOOP - Structural Analysis with Object-Oriented Programming
% This is the main script file of StAnOOP.
% Run this script to select the input file with model and analysis data for
% performing the simulation.
% Multiple files can be selected, by holding down the Shift or Ctrl key and
% clicking file names, to run the simulations sequentially.

%% Plotting Options (NOT WORKING YET)
% Set plotting options with flags 0 (NO) or 1 (YES):
%   eid -> Plot element numbers
%   nid -> Plot node numbers
%   scl -> Scale factor for deformed mesh
%   dx  -> Plot countours of displacements in X direction
%   dy  -> Plot countours of displacements in Y direction
%   sx  -> Plot countours of normal stresses in X direction
%   sy  -> Plot countours of normal stresses in Y direction
%   txy -> Plot countours of shear stresses
%   s1  -> Plot countours of major principal stresses
%   s2  -> Plot countours of minor principal stresses
%   tmx -> Plot countours of maximum shear stresses

% Mesh options
eid = 1;
nid = 1;

% Result options:
scl = 1.0;
dx  = 1;   dy  = 1;
sx  = 1;   sy  = 1;   txy = 1;
s1  = 1;   s2  = 1;   tmx = 1;

%% Run Simulation
opt = [eid,nid,scl,dx,dy,sx,sy,txy,s1,s2,tmx];
clc; close all; clearvars -except opt; 
addpath('src');
sim = Simulation();
sim.runAll(opt);