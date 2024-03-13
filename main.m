%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elasticity2D - FEM                                                    %
% Main driver file.                                                     % 
% This file contains MATLAB code of a program for linear-elastic,       %
% displacement-based, two-dimensional, finite-element analysis for      %
% solving a stress-distribution elasticity problem.                     % 
% The program reads a file with FE model data, in a neutral format,     %
% assembles a system of equations, solves the system and visualizes     %
% the response data.                                                    %
%                                                                       %
% Author:                                                               %
% Luiz Fernando Martha                                                  %
% Pontifical Catholic University of Rio de Janeiro - PUC-Rio            %
% Department of Civil Engineering and Tecgraf                           %
%                                                                       %
% Adapted from the program Heat2D                                       %
% by Haim Waisman, Rensselaer Polytechnic Institute.                    %
% Available in the the companion site (http://1coursefem.blogspot.com)  %
% of the book by Fish, J. and Belytschko, T., A First Course in Finite  %
% Elements - Chapter 12: Finite Element Programming with MATLAB, 2007.  %
%                                                                       %
% It is assumed that there is only one type of element and only one     %
% type of gauss integration order in each finite element model.         %
%                                                                       %
% It is assumed that there is a single load case.                       %
%                                                                       %
% See global variables in file include_gbl_refs.m.                      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all; 

% Include global variables
include_gblrefs;  

% Initialize constant global variables
init_constants;  

% Preprocessing
[K,F,D] = preProcessor;

% Assemble global stiffness matrix
fprintf(1,'Assembling stiffness matrix...\n');
for e = 1:nel
 ke = elemStiffMtx(e);
 K = drvAssembleMtx(K,e,ke);
end

% Assemble global forcing vector
fprintf(1,'Assembling forcing vector...\n');
F = drvPointLoads(F);    % initialize forcing vector with nodal point loads
F = drvEdgeLoads(F);     % add edge (element side) loads to forcing vector
F = drvAreaLoads(F);     % add area (element) loads to forcing vector

% Partition and solve system of equations
fprintf(1,'Solving system of equations...\n');
[D,F] = drvSolveEqnSystem(neq,neqfree,K,F,D);

% Postprocessing
posProcessor(D);

fprintf(1,'Finished.\n');
