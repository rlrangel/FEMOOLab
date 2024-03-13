# FEMOOLab: Finite Element Method Object-Oriented Laboratory

<p align=center><img height="100.0%" width="100.0%" src="https://github.com/rlrangel/FEMOOLab/blob/master/docs/images/logo.png"></p>

[![Release][release-image]][release] [![License][license-image]][license] [![Contributing][contributing-image]][contributing]

[release-image]: https://img.shields.io/badge/release-1.0.0-green.svg?style=flat
[release]: https://github.com/rlrangel/FEMOOLab/releases

[license-image]: https://img.shields.io/badge/license-MIT-green.svg?style=flat
[license]: https://github.com/rlrangel/FEMOOLab/blob/master/LICENSE

[contributing-image]: https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg
[contributing]: https://github.com/rlrangel/FEMOOLab/blob/master/CONTRIBUTING.md

FEMOOLab is a MATLAB program for performing FEM-based numerical simulations, implemented in a modular OOP framework to allow different types of models and physics.

## Table of Contents
- [Main Features](#main-features)
- [Implementation Aspects](#implementation-aspects)
- [Instructions](#instructions)
    - [Input Files](#input-files)
    - [Running Simulations](#running-simulations)
    - [Testing](#testing)
- [Examples](#examples)
- [Documentation](#documentation)
- [How to Contribute](#how-to-contribute)
- [How to Cite](#how-to-cite)
- [Authorship](#authorship)
- [Acknowledgement](#acknowledgement)
- [License](#license)

## Main Features

The program solves 2D steady-state and transient problems of
structural analysis (linear-elasticity) and
thermal analysis (conductive and convective heat transfer)
with isoparametric finite element formulations.

The available **analysis model types** are:
- Structural plane stress
- Structural plane strain
- Structural axisymmetric
- Structural thick (mindlin) plate
- Thermal plane conduction
- Thermal axisymmetric conduction
- Thermal plane convection-diffusion

The available **element types** are:
- T3: Linear planar triangular isoparametric element with 3 nodes and Lagrangean interpolation
- Q4: Linear planar quadrilateral isoparametric element with 4 nodes and Lagrangean interpolation
- T6: Quadratic planar triangular isoparametric element with 6 nodes and Lagrangean interpolation
- Q8: Quadratic planar quadrilateral isoparametric element with 8 nodes and Serendipity interpolation

## Implementation Aspects

FEMOOLab is fully written in the [MATLAB][matlab_website] programming language,
and adopts the Object Oriented Programming (OOP) paradigm to offer modularity and extensibility.

The source code can run in any operating system where MATLAB can be installed.

## Instructions

### Input Files

The program reads a file with finite element model data that follows Tecgraf's [neutral file][nf_link] format.

### Running Simulations

To run a simulation, launch MATLAB and execute the script file [*main.m*][main_file_link] located inside the folder [*src*][src_folder_link].

A dialog box will pop up to select appropriate input files.
Multiple input files can be selected to run simulations sequentially, as long as they are located in the same directory.

If the models and parameters are read correctly, the simulations are started and their progresses are printed in the MATLAB command window.

Graphical results are plotted in MATLAB figure windows.

### Testing

In progress...

## Examples

A collection of sample models are available inside the folder [*examples*][examples_link].

## Documentation

In progress...

## How to Contribute

Please check the [contribution guidelines][contribute_link].

## How to Cite

To cite this repository, you can use the metadata from [this file][citation_link].

## Authorship

- **Rafael Rangel**<sup>1</sup> (<rrangel@cimne.upc.edu>)
- **Luiz Fernando Martha**<sup>2,3</sup> (<lfm@tecgraf.puc-rio.br>)
- **João Carlos Peixoto**<sup>2,3</sup> (<aaaaaaaa@tecgraf.puc-rio.br>)

<sup>1</sup> International Center for Numerical Methods in Engineering ([CIMNE][cimne_website])

<sup>2</sup> Pontifical Catholic University of Rio de Janeiro (PUC-Rio) - [Department of Civil and Environmental Engineering][civil_website]

<sup>3</sup> Tecgraf Institute of Technical-Scientific Software Development of PUC-Rio ([Tecgraf/PUC-Rio][tecgraf_website])

## License

FEMOOLab is licensed under the [MIT license][mit_license_link],
which allows the program to be freely used by anyone for modification, private use, commercial use, and distribution, only requiring preservation of copyright and license notices.
No liability and warranty are provided.

[matlab_website]:   https://www.mathworks.com/
[nf_link]:          https://web.tecgraf.puc-rio.br/neutralfile
[main_file_link]:   https://github.com/rlrangel/FEMOOLab/blob/master/src/main.m
[src_folder_link]:  https://github.com/rlrangel/FEMOOLab/tree/master/src
[examples_link]:    https://github.com/rlrangel/FEMOOLab/tree/master/examples
[contribute_link]:  https://github.com/rlrangel/FEMOOLab/blob/master/CONTRIBUTING.md
[citation_link]:    https://github.com/rlrangel/FEMOOLab/blob/master/CITATION.cff
[mit_license_link]: https://choosealicense.com/licenses/mit/
[cimne_website]:    https://www.cimne.com/
[civil_website]:    https://www.civ.puc-rio.br/en/
[tecgraf_website]:  https://www.tecgraf.puc-rio.br/
