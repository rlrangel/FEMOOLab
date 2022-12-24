%% FEMOOLab - Finite Element Method Object-Oriented Laboratory
%
%% Instructions
%
% This is the main script file of the FEMOOLab program.
%
% To run a simulation, execute this script and select an appropriate
% parameters file with the _.json_ extension.
%
% Multiple parameter files can be selected to run simulations sequentially,
% as long as they are located in the same folder.
%
% Sub-folders with the simulation name plus the suffix "_out" are created
% to receive the output files with the results of each simulation.
%
%% OOP classes
%
% This program adopts an Object Oriented Programming (OOP) paradigm.
% The following OOP super-classes are implemented:
%
% TODO...
%
%% Initialization
clc; clearvars; close all;
addpath(genpath(pwd));
Master().RunSimulations(1);