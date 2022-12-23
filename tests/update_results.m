%% Regression Tests
%
% This script file updates the reference result files of the selected test
% models with the current obtained results.
%
clc; clearvars; close all; format long;

% Decimal precision for printing current results into reference files
% (e.g. 1e-prec)
prec = 16;

% Run tests
addpath(genpath(pwd));
addpath(genpath('..\src'));
run(2,prec);