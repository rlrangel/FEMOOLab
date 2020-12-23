%% StAnOOP - Finite Element Model Laboratory
%
%% Description
%
% This is the main script file of StAnOOP. This is a MATLAB program for
% linear-elastic, displacement-based, static analysis of two-dimensional
% finite elements structure models.
%
% For each structural analysis, the program assembles a system of
% equilibrium equations, solves the system and displays the analysis
% results.
%
%% Usage
% Run this script to select the input file with model and analysis data for
% performing the simulation.
% Multiple files can be selected, by holding down the Shift or Ctrl key and
% clicking file names, to run the simulations sequentially.
%
%% Brief description
% This version (1.1) of the program handles plane models with C0 continuity
% for conventional isoparametric finite element types. 
%
% The Object-Oriented Programming (OOP) paradigm is adopted in the
% implementation of the analysis process. The use of OOP is justified by
% the clarity, organization and generality that this programming paradigm
% provides to the code.
%
% This is a non-graphical version of program that may be used in a extended
% to a GUI (Graphical User Interface) version.
% The non-graphical version reads a structural model from a neutral
% format file and prints displays analysis results in separate MATLAB
% figures.
%
%% Authors
%%%
% * Rafael Lopez Rangel (rafaelrangel@tecgraf.puc-rio.br)
% * Luiz Fernando Martha (lfm@tecgraf.puc-rio.br)
%
% Pontifical Catholic University of Rio de Janeiro - PUC-Rio
%
% Department of Civil and Environmental Engineering and Tecgraf Institute
% of Technical-Scientific Software Development of PUC-Rio (Tecgraf/PUC-Rio)
%
%% History
% @version 1.1
%
%%%
% * Initial version 1.0: August 2018
%
% Initially prepared by Rafael Lopez Rangel as a semester project for the
% course CIV 2118 - Método dos Elementos Finitos, second term, Department
% of Civil and Environmental Engineering, PUC-Rio.
%
% Based on the program Elasticity2D (non OOP) by Luiz Fernando Martha and 
% on LESM version 1, developed in the undergraduate thesis "Development of
% a Graphic Program for Structural Analysis of Linear Element Models", by
% Rafael Lopez Rangel, advised by Professor Luiz Fernando Martha,
% Department of Civil and Environmental Engineering, PUC-Rio, December, 2016.
%
%%%
% * Version 1.1: October 2020
%
% Version 1.1 was prepared by Luiz Fernando Martha for the course CIV 2801 -
% Fundamentos de Computação Gráfica Aplicada.
%
% In the version, Element class was split into three classes: Element,
% Shape, and IntegPt.
%
% Input data files now follow specification of Tecgraf's neutral file 
% format: <https://web.tecgraf.puc-rio.br/neutralfile NF spec>
%
%% Analysis Model Types
%
% Current version of StAnOOP program considers only static linear-elastic
% structural analysis of 2D (plane) finite elements models of C0 continutity,
% which can be of the following types:
%%%
% * Plane stress analysis
% * Plane strain analysis
% * Axi-symmetric analysis
%
%% Element Shape Types
%
% Current version of StAnOOP adopts an isoparametric finite element
% formulation with conventional element shape types:
%%%
% * Tria3: linear planar triangular isoparametric element with 3 nodes,
%          with Lagrangean interpolation
% * Quad4: linear planar quadrilateral isoparametric element with 4 nodes,
%          with Lagrangean interpolation
% * Tria6: quadratic planar triangular isoparametric element with 6 nodes,
%          with Lagrangean interpolation
% * Quad8: quadratic planar quadrilateral isoparametric element with 8 nodes,
%          with Serendipity interpolation
%
%% Load Type
%
% There are four types of loads considered in the StAnOOP program:
%%%
% * Concentrated nodal force and moment in global axes directions.
% * Uniformely distributed force along a side of a planar or solid element.
% * Uniformely distributed force along a planar element entire domain.
%
% In addition, nodal prescribed displacements and rotations may be specified.
%
% It is assumed that there is a single load case.
%
%% Materials
%
% All materials in StAnOOP are considered to have linear elastic behavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in
% all directions.
%
%% Thicknesses
%
% A list of thickness properties for planar element is provided.
%
%% Object Oriented Classes
%
% This program adopts an Object Oriented Programming (OOP) paradigm, in
% which the following OOP classes are implemented:
%%%
% * <simulation.html Simulation: simulation driver class>
% * <read.html Read: Read input data driver class>
% * <result.html Result: Plot results driver class>
% * <anl.html Anl: analysis class>
% * <model.html Model: model class>
% * <anm.html Anm: analysis model class>
% * <node.html Node: node class>
% * <element.html Element: finite element class>
% * <shape.html Shape: element shape class>
% * <material.html Material: material class>
% * <gauss.html Gauss: integration quadrature class>
%
