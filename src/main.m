%% FEMOOLab - Finite Element Method Object-Oriented Laboratory
%
% <<logo_femoolab.png>>
%
% This program uses the Finite Element Method (FEM) to simulate
% multiphysics problems in steady-state or transient regime.
%
% Currently, linear analyses of structural and thermal models are available.
%
% The program is designed to be an educational tool, aimed at those who are
% learning the implementation aspects about the FEM.
%
% For more information, visit the
% <https://gitlab.com/rafaelrangel/femoolab GitLab repository>.
%
%% Instructions
%
% This is the main script file of the FEMOOLab program.
%
% To *run a simulation*, execute this script and select an appropriate
% parameters file with the _.json_ extension.
%
% To *load results* from a previously run simulation, execute this
% script and select an appropriate storage file with the _.mat_ extension.
%
% There are two ways to select an input file:
%
% * Provide its complete path thorugh the string variable _file_ in the end
%   of this script (e.g. file = 'C:\...\ProjectParameters.json').
% * Leave the string variable _file_ blank (file = '') and use the
%   selection dialog that appears when running this script (multiple files
%   can be selected to run simulations sequentially).
%
% A folder with the problem name plus the suffix "_out" is created to
% receive the output files.
%
%% OOP classes
%
% This program adopts an Object Oriented Programming (OOP) paradigm.
% The following OOP super-classes are implemented:
%
% * <master.html Master>
%
%% Authors
%
% * Rafael Rangel (rrangel@cimne.upc.edu)
%
% * Luiz Fernando Martha (lfm@tecgraf.puc-rio.br)
%
% International Center for Numerical Methods in Engineering
% (<https://www.cimne.com/ CIMNE>)
% and 
% Tecgraf Institute of Technical-Scientific Software Development of PUC-Rio
% (<https://www.tecgraf.puc-rio.br/ Tecgraf/PUC-Rio>)
%
% <<logo_institutions.png>>
%
%% History
%
% * Version 1.0 - December, 2021
%
%% Initialization
clc; clearvars; close all;
file = '';
addpath('gauss',...
        'shape',...
        'utl');
Master().execute(file);