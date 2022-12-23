# FEMOOLab - Finite Element Method Object-Oriented Laboratory

<p align=center><img height="100.0%" width="100.0%" src="https://gitlab.com/rafaelrangel/femoolab/-/raw/new-version/docs/images/logos/logo_femoolab.png"></p>

FEMOOLab is a MATLAB program for performing numerical simulations using the Finite Element Method (FEM).

Its purpose is to offer an open-source tool that is modular and extensible, allowing easy implementation, inspection and testing of new models and formulations.
It targets researchers, who can use it as a numerical laboratory, and teachers and students who want to teach or learn implementation aspects about the FEM.

The program is based on the [FEMOOP][femoop_website] system.

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
- [License](#license)

## Main Features

**TODO**

The program performs thermomechanical analyses of two-dimensional (2D) models.

Steady-state and transient analysis.

Linear-elastic analysis.

The program handles plane models with C0 continuity for conventional isoparametric finite element types. 

Analysis Model Types: 2D (plane) finite elements models of C0 continutity.

The program adopts an isoparametric finite element formulation with conventional element shape types:
* _Tria3_: linear planar triangular isoparametric element with 3 nodes, with Lagrangean interpolation
* _Quad4_: linear planar quadrilateral isoparametric element with 4 nodes, with Lagrangean interpolation
* _Tria6_: quadratic planar triangular isoparametric element with 6 nodes, with Lagrangean interpolation
* _Quad8_: quadratic planar quadrilateral isoparametric element with 8 nodes, with Serendipity interpolation

## Implementation Aspects

FEMOOLab is fully written in the [MATLAB][matlab_website] programming language,
and adopts the Object Oriented Programming (OOP) paradigm to offer modularity and extensibility.

The source code can run in any operating system where MATLAB can be installed
(the program is tested for version 2022a of MATLAB).

Because it is developed with a high-level interpreted programming language using serial processing,
code efficiency is not a priority and therefore only small to medium-scale problems should be simulated.

## Instructions

### Input Files

Two files are used as input to run a simulation:

* **Model Parts (_.txt_)**: 

Text file with node IDs and coordinates, and element IDs and connectivity.
It also groups nodes and elements into model parts, which are used to apply attributes and conditions to multiple nodes and/or elements from the _Project_ _Parameters_ file.

The file name must be indicated in the accompanying _Project_ _Parameters_ file.

A tutorial on this file can be found on its [Wiki page][wiki_mparts_link].
Moreover, a [template][modelparts_link] of this file, with all the possible input options, is available.

* **Project Parameters (_.json_)**: 

[JSON][json_link] file with all the parameters and options for the analysis and outputs, as well as the attributes and conditions applied to the model.

The name of the accompanying _Model_ _Parts_ file must be indicated in the Input Group ["ProblemData"][problem_data_link] of this file.

A tutorial explaining each input field of this file can be found on its [Wiki page][wiki_parameters_link].
Moreover, a [template][modelparts_link] of this file, with all the possible input options, is available.

### Running Simulations

To run a simulation, launch MATLAB and execute the script file [*main.m*][main_file_link] located inside the folder [*src*][src_folder_link].

A dialog box will pop up to select an appropriate _Project_ _Parameters_ file.
Multiple _Project_ _Parameters_ files can be selected to run simulations sequentially, as long as they are located in the same directory.

If the models and parameters are read correctly, the simulations are started and their progresses are printed in the MATLAB command window.

Sub-folders with the names of the simulations, plus the suffix _"out"_, are created to receive the output files with the requested results of each simulation.

### Testing

Recursive tests are available to verify that the program is working correctly and that the current results are matching with the reference results.
The reference results are stored in files with a _.pos_ extension.

To **run the tests**, execute the script file [*run_tests.m*][run_tests_link] located inside the folder [*tests*][tests_folder_link].
In that file, you can set the tolerance for comparing the current results with the reference results.
A dialog box will pop up to select the _Project_ _Parameters_ files of the tests to be run, located inside the sub-folder [*test_models*][test_models_link].
The result of each test is then printed in the MATLAB command window.

To **generate or update the reference results** of the test models, execute the script file [*update_results.m*][upd_tests_link] located inside the folder [*tests*][tests_folder_link].
In that file, you can set the number of decimal places to print the reference results.
A dialog box will pop up to select the _Project_ _Parameters_ files of the tests to be updated, which are located inside the sub-folder [*test_models*][test_models_link].
It will run the selected tests and overwrite existing reference results.

## Examples

A collection of sample models is available inside the folder [*examples*][examples_link].

The examples are separated into different sub-folders according to their analysis type,
and each example has its _Project_ _Parameters_ and _Model_ _Parts_ files, as well as some results in the output sub-folders.

## Documentation

General and important information for starting using the program is available in the [Wiki pages][wiki_link].

The OOP classes of the program are documented in [UML diagrams][uml_link] and their codes in HTML files, located inside the folder [*html*][html_folder_link].
These files can be browsed on their [Wiki page][wiki_html_link].

Tutorials explaining all the options of the program that can be added to the input files can also be found on the [Wiki page][wiki_link].

## How to Contribute

Please check the [contribution guidelines][contribute_link].

## How to Cite

**TODO**

## Authorship

- **Luiz Fernando Martha** (<lfm@tecgraf.puc-rio.br>)
- **Rafael Rangel** (<rrangel@cimne.upc.edu>)
- **João Carlos Peixoto** (<joaoc_peixoto@hotmail.com>)

Pontifical Catholic University of Rio de Janeiro (PUC-Rio) -
[Department of Civil and Environmental Engineering][civil_website]

Tecgraf Institute of Technical-Scientific Software Development of PUC-Rio
([Tecgraf/PUC-Rio][tecgraf_website])

International Center for Numerical Methods in Engineering ([CIMNE][cimne_website]) 

<p float="left">
&nbsp;&nbsp;&nbsp;&nbsp;
<img src="https://gitlab.com/rafaelrangel/femoolab/-/raw/new-version/docs/images/logos/logo_puc.png" width="200"/>
&nbsp;&nbsp;&nbsp;&nbsp;
<img src="https://gitlab.com/rafaelrangel/femoolab/-/raw/new-version/docs/images/logos/logo_tecgraf.png" width="280"/> 
&nbsp;&nbsp;&nbsp;&nbsp;
<img src="https://gitlab.com/rafaelrangel/femoolab/-/raw/new-version/docs/images/logos/logo_cimne.png" width="240"/>
</p>

## License

FEMOOLab is licensed under the [MIT license][mit_license_link],
which allows the program to be freely used by anyone for modification, private use, commercial use, and distribution, only requiring preservation of copyright and license notices.

No liability and warranty is provided.
Neither the authors nor any related institution are responsible for any use or misuse of the program and the results.
The aforementioned assume no liability or responsibility to any person or company for direct or indirect damages resulting from the use of any information or the use of any of the software made available here.
The user is responsible for any and all conclusions made while using the program.

[femoop_website]:       http://webserver2.tecgraf.puc-rio.br/femoop/
[matlab_website]:       https://www.mathworks.com/
[wiki_mparts_link]:     https://gitlab.com/rafaelrangel/femoolab/-/wikis/model-parts-file
[modelparts_link]:      https://gitlab.com/rafaelrangel/femoolab/-/blob/master/docs/help/ModelParts_template.txt
[parameters_link]:      https://gitlab.com/rafaelrangel/femoolab/-/blob/master/docs/help/ProjectParameters_template.json
[json_link]:            https://www.json.org/
[problem_data_link]:    https://gitlab.com/rafaelrangel/femoolab/-/wikis/project-parameters-file#problemdata
[wiki_parameters_link]: https://gitlab.com/rafaelrangel/femoolab/-/wikis/project-parameters-file
[main_file_link]:       https://gitlab.com/rafaelrangel/femoolab/-/blob/master/src/main.m
[src_folder_link]:      https://gitlab.com/rafaelrangel/femoolab/-/tree/master/src
[run_tests_link]:       https://gitlab.com/rafaelrangel/femoolab/-/blob/master/tests/run_tests.m
[tests_folder_link]:    https://gitlab.com/rafaelrangel/femoolab/-/tree/master/tests
[test_models_link]:     https://gitlab.com/rafaelrangel/femoolab/-/tree/master/tests/test_models
[upd_tests_link]:       https://gitlab.com/rafaelrangel/femoolab/-/blob/master/tests/update_results.m
[examples_link]:        https://gitlab.com/rafaelrangel/femoolab/-/tree/master/examples
[wiki_link]:            https://gitlab.com/rafaelrangel/femoolab/-/wikis/home
[uml_link]:             https://www.uml.org/
[html_folder_link]:     https://gitlab.com/rafaelrangel/femoolab/-/tree/master/docs/html
[wiki_html_link]:       https://gitlab.com/rafaelrangel/femoolab/-/wikis/classes-documentation
[contribute_link]:      https://gitlab.com/rafaelrangel/femoolab/-/blob/master/CONTRIBUTING.md
[civil_website]:        https://bananastudio.com.br/civil/web/
[tecgraf_website]:      https://www.tecgraf.puc-rio.br/
[cimne_website]:        https://www.cimne.com/
[mit_license_link]:     https://choosealicense.com/licenses/mit/