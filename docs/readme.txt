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
% It is assumed that all nodes and elements are numbered consecutively from
% one to the total number of nodes and elements.
% It is assumed that there is only one type of element and only one
% type of gauss integration order in each finite element model.
% It is assumed that there is a single load case.

===========================================================================
%% StAnOOP - Structural Analysis with Object-Oriented Programming
%
% This is the main file of *StAnOOP*. This is a MATLAB program for
% linear-elastic, displacement-based, static analysis of bi-dimensional
% and tri-dimensional linear element models, and bi-dimensional finite
% element models for solving a stress-distribution elasticity problem.
% For each structural analysis, the program uses the direct stiffness method
% to assemble and solve the system of equilibrium equations.
% The program reads a structural model data from a neutral format
% file and prints model information and analysis results in a text file.
% The program adopts an object-oriented programming (OOP) paradigm, so
% it is composed by classes responsible for instantiating objects to
% perform the structural analysis.
%
%% Authors
%%%
% * Luiz Fernando Martha (lfm@tecgraf.puc-rio.br)
%%%
% * Rafael Lopez Rangel (rafaelrangel@tecgraf.puc-rio.br)
%%%
% Pontifical Catholic University of Rio de Janeiro - PUC-Rio
%
% Department of Civil and Environmental Engineering and Tecgraf Institute
% of Technical-Scientific Software Development of PUC-Rio (Tecgraf/PUC-Rio)
%
%% Background
%
% Initially developed as an assignment for the course CIV 2118 (Introducao
% ao Metodo dos Elementos Finitos), 2016, second term, by professor
% Deane Roehl - Department of Civil and Environmental Engineering, PUC-Rio.
%
% Adapted from:
%%%
% * *Heat2D*: by Haim Waisman, Rensselaer Polytechnic Institute. 
% Available in the the companion website http://1coursefem.blogspot.com
% of the book by Fish, J. and Belytschko, T., A First Course in Finite
% Elements - Chapter 12: Finite Element Programming with MATLAB, 2007.
%%%
% * *LESM*: by Rafael Rangel and Luiz Fernando Martha, PUC-Rio. 
% Available in the the companion website <https://web.tecgraf.puc-rio.br/lesm/>
%
%% Model types
%
% * *Linear Element Models*:
% Linear element models are made of uniaxial elements (bars and beams),
% with one dimension much larger than the others, such as trusses and frames.
% The analysis process of this type of model is based on a nodal
% dicretization, where a global analysis result is computed based on nodal
% displacements and rotations and a local analysis result is computed with
% the effect of loads inside each element. These analyzes uses analytical
% expressions of linear element behavior. Therefore, the analytical
% solution of the problem is then obtained by the sum of the two analyzes.
%
% * *Finite Element Models*:
% Bi-dimensional finite element models are made of plane elements with one
% dimension much smaller than the others, such as plates.
% The analysis process of this type of model considers only the effect of
% nodal displacements, neglecting local effects. The solution of nodal
% displacements are based on approximate expressions for the element
% behavior. Therefore, the solution of the problem is approximated.
%
%% Analysis model types
%
% *StAnOOP* considers only static linear-elastic analysis of
% linear element models and bi-dimensional finite element models.
%%%
% *Linear Element Models*
%
% * Plane Truss
% * Plane Frame
% * Grillage
% * Spatial Truss
% * Spatial Frame
%
% For linear element analysis models, the program reads a neutral format
% file with the _.lsm_ extension, produced by the *LESM* program.
%
%%%
% *Finite Element Models*
%
% * Plane Stress
% * Plane Strain
% * Axisymmetric
%
% For finite element analysis models, the program reads a neutral format
% file with the _.nf_,_.pos_ or _.dat_ extensions, produced by the
% *Sigma2D* program.
%
%% Element types
%
% *StAnOOP* assumes that there is only one type of element in each model and
% only one type of gauss integration order in finite element models.
%%%
% *Linear Element Models*
%
% * Navier (Euler-Bernoulli) beam element.
% * Timoshenko beam element
%
%%%
% *Finite Element Models*
%
% * 3 node triangular element (T3)
% * 4 node rectangular element (Q4)
% * 6 node triangular element (T6)
% * 8 node rectangular element (Q8)
%
%% Load types
%
% It is assumed that there is a single load case.
%%%
% *Linear Element Models*
%
% * Concentrated nodal load in global axes directions.
% * Uniformely distributed force on elements, spanning its entire
%   length, in local or in global axes directions.
% * Linearly distributed force on elements, spanning its entire length,
%   in local or in global axes directions.
% * Uniform temperature variation on faces of elements.
%
% In addition, nodal prescribed displacement and rotations may be
% specified.
%
%%%
% *Finite Element Models*
%
% * Concentrated nodal load.
% * Uniformely distributed edge load on elements.
% * Uniformely distributed area load on elements.
%
% All loads in finite element models are given in global axes directions.
% In addition, nodal prescribed displacement may be specified.
%
%% Geometry
%
% In linear element models, all cross-sections are considered to be uniform
% and of a generic type, which means that their shapes are not specified,
% only their geometric properties are provided, such as area, moment of
% inertia and height.
%
% In finite element models, all thicknesses are constant with strain,
% stresses, and internal forces varying uniformly along the thickness.
%
%% Materials
%
% All materials in *StAnOOP* are considered to have linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in
% all directions.
%
%% Object oriented classes
% This program adopts an Object Oriented Programming (OOP) paradigm, in
% which the following OOP classes are implemented:
%%%
% * <model.html Model: model class>.
% * <anm.html Anm: analysis model class>.
% * <node.html Node: node class>.
% * <elem.html Elem: element class>.
% * <lelem.html Lelem: load element class>.
% * <material.html Material: material class>.
% * <section.html Section: cross-section class>.
% * <print.html Print: print class>.
%
%% Auxiliary functions and files
%%%
% * <include_constants.html include_constants: File with global constants>.
% * <preReadLSM.html preReadLSM: Reads linear elements structural model data from a neutral format file>.
% * <preReadNF.html preReadNF: Reads finite element model data from a neutral format file>.
%