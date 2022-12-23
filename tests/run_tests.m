%% Regression Tests
%
% This script file runs the selected test models and checks if the obtained
% results agree with the reference results.
%
clc; clearvars; close all; format long;

% Tolerance for comparing current results with reference
tol = 1e-12;

% Run tests
addpath(genpath(pwd));
addpath(genpath('..\src'));
run(1,tol);