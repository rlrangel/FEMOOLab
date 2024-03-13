# FEMOOLab: Finite Element Method Object-Oriented Laboratory

Finite Element Method Object-Oriented Laboratory (FEMOOLab)
is a MATLAB program to perform FEM-based numerical simulations,
implemented in a modular OOP framework to allow different types of models and physics.

## Table of Contents


## Main Features

linear-elastic, displacement-based, static analysis of two-dimensional
finite elements structure models.

% This is the main file of *StAnOOP*. This is a MATLAB program for
% linear-elastic, displacement-based, static analysis of bi-dimensional
% and tri-dimensional linear element models, and bi-dimensional finite
% element models for solving a stress-distribution elasticity problem.

% This file contains MATLAB code of a program for linear-elastic,       %
% displacement-based, two-dimensional, finite-element analysis for      %
% solving a stress-distribution elasticity problem.      

This version (1.1) of the program handles plane models with C0 continuity
for conventional isoparametric finite element types. 

This is a non-graphical version of program that may be used in a extended
to a GUI (Graphical User Interface) version.
The non-graphical version reads a structural model from a neutral
format file and prints displays analysis results in separate MATLAB
figures.

Analysis Model Types

Current version of FEMOOLab program considers only static linear-elastic
structural analysis of 2D (plane) finite elements models of C0 continutity,
which can be of the following types:
* Plane stress analysis
* Plane strain analysis
* Axi-symmetric analysis

%% Analysis model types
% *StAnOOP* considers only static linear-elastic analysis of
% linear element models and bi-dimensional finite element models.

% *Finite Element Models*
%
% * Plane Stress
% * Plane Strain
% * Axisymmetric

%% Element types
%
% *StAnOOP* assumes that there is only one type of element in each model and
% only one type of gauss integration order in finite element models.

% *Finite Element Models*
%
% * 3 node triangular element (T3)
% * 4 node rectangular element (Q4)
% * 6 node triangular element (T6)
% * 8 node rectangular element (Q8)

%% Load types
% It is assumed that there is a single load case.
% * Concentrated nodal load.
% * Uniformely distributed edge load on elements.
% * Uniformely distributed area load on elements.
% All loads in finite element models are given in global axes directions.
% In addition, nodal prescribed displacement may be specified.

Element Shape Types

Current version of FEMOOLab adopts an isoparametric finite element
formulation with conventional element shape types:
* Tria3: linear planar triangular isoparametric element with 3 nodes,
         with Lagrangean interpolation
* Quad4: linear planar quadrilateral isoparametric element with 4 nodes,
         with Lagrangean interpolation
* Tria6: quadratic planar triangular isoparametric element with 6 nodes,
         with Lagrangean interpolation
* Quad8: quadratic planar quadrilateral isoparametric element with 8 nodes,
         with Serendipity interpolation

Load Type

There are four types of loads considered in the FEMOOLab program:
* Concentrated nodal force and moment in global axes directions.
* Uniformely distributed force along a side of a planar or solid element.
* Uniformely distributed force along a planar element entire domain.

In addition, nodal prescribed displacements and rotations may be specified.
It is assumed that there is a single load case.

%% Geometry
% In finite element models, all thicknesses are constant with strain,
% stresses, and internal forces varying uniformly along the thickness.

Materials

All materials in FEMOOLab are considered to have linear elastic behavior.
In adition, homogeneous and isotropic properties are also considered,
that is, all materials have the same properties at every point and in
all directions.

% All materials in *StAnOOP* are considered to have linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in
% all directions.

Thicknesses

A list of thickness properties for planar element is provided.

## Implementation Aspects

The Object-Oriented Programming (OOP) paradigm is adopted in the
implementation of the analysis process. The use of OOP is justified by
the clarity, organization and generality that this programming paradigm
provides to the code.

% The program adopts an object-oriented programming (OOP) paradigm, so
% it is composed by classes responsible for instantiating objects to
% perform the structural analysis.

Object Oriented Classes

This program adopts an Object Oriented Programming (OOP) paradigm, in
which the following OOP classes are implemented:
* <simulation.html Simulation: simulation driver class>
* <read.html Read: Read input data driver class>
* <result.html Result: Plot results driver class>
* <anl.html Anl: analysis class>
* <model.html Model: model class>
* <anm.html Anm: analysis model class>
* <node.html Node: node class>
* <element.html Element: finite element class>
* <shape.html Shape: element shape class>
* <material.html Material: material class>
* <gauss.html Gauss: integration quadrature class>

## Instructions

### Input Files

Input data files now follow specification of Tecgraf's neutral file 
format: <https://web.tecgraf.puc-rio.br/neutralfile NF spec>

% The program reads a file with FE model data, in a neutral format, 

% It is assumed that there is only one type of element and only one     %
% type of gauss integration order in each finite element model.         %
%                                                                       %
% It is assumed that there is a single load case.                       %

% It is assumed that all nodes and elements are numbered consecutively from
% one to the total number of nodes and elements.
% It is assumed that there is only one type of element and only one
% type of gauss integration order in each finite element model.
% It is assumed that there is a single load case.

### Running Simulations

Run this script to select the input file with model and analysis data for
performing the simulation.
Multiple files can be selected, by holding down the Shift or Ctrl key and
clicking file names, to run the simulations sequentially.

### Testing

## Examples

## Documentation



## How to Contribute



## How to Cite



## Authorship

* Rafael Lopez Rangel (rafaelrangel@tecgraf.puc-rio.br)
* Luiz Fernando Martha (lfm@tecgraf.puc-rio.br)

Pontifical Catholic University of Rio de Janeiro - PUC-Rio

Department of Civil and Environmental Engineering and Tecgraf Institute
of Technical-Scientific Software Development of PUC-Rio (Tecgraf/PUC-Rio)

## Acknowledgement



## License

