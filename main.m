%% FEMOOLab - Finite Element Method Object-Oriented Laboratory
% This is the main script file of FEMOOLab.
% Run this script to select the input file with model and analysis data for
% performing the simulation.
% Multiple files can be selected, by holding down the Shift or Ctrl key and
% clicking file names, to run the simulations sequentially.
% For more information, see the README file in the documents folder.

%% Plotting Options
% Set result plotting options with flags (true or false):
% * eid    -> Plot element numbers (NOT IMPLEMENTED)
% * nid    -> Plot node numbers (NOT IMPLEMENTED)
% * scl    -> Scale factor for deformed mesh (NOT IMPLEMENTED)
% * dx     -> Plot contour of displacements in X direction (NOT IMPLEMENTED)
% * dy     -> Plot contour of displacements in Y direction (NOT IMPLEMENTED)
% * dz     -> Plot contour of displacements in Z directiont (NOT IMPLEMENTED)
% * rx     -> Plot contour of rotations about X axis (NOT IMPLEMENTED)
% * ry     -> Plot contour of rotations about Y axis (NOT IMPLEMENTED)
% * rz     -> Plot contour of rotations about Z axis (NOT IMPLEMENTED)
% * smooth -> Smooth element results at common nodes
% * sxx    -> Plot contour of normal stresses in X direction
% * syy    -> Plot contour of normal stresses in Y direction
% * szz    -> Plot contour of normal stresses in Z direction
% * txy    -> Plot contour of XY shear stresses
% * txz    -> Plot contour of XZ shear stresses
% * tyz    -> Plot contour of YZ shear stresses
% * s1     -> Plot contour of principal stresses 1
% * s2     -> Plot contour of principal stresses 2
% * s3     -> Plot contour of principal stresses 3
% * taumax -> Plot contour of maximum shear stresses
% * mxx    -> Plot contour of moment about X direction (NOT IMPLEMENTED)
% * myy    -> Plot contour of moment about Y direction (NOT IMPLEMENTED)
% * mxy    -> Plot contour of torsion moment (NOT IMPLEMENTED)
% * qxz    -> Plot contour of XZ shear force (NOT IMPLEMENTED)
% * qyz    -> Plot contour of YZ shear force (NOT IMPLEMENTED)
% * m1     -> Plot contour of principal moment 1 (NOT IMPLEMENTED)
% * m2     -> Plot contour of principal moment 2 (NOT IMPLEMENTED)
% * tormax -> Plot contour of maximum torsion (NOT IMPLEMENTED)

%% Initialization

% Clear workspace
clear
close(findall(0,'Type','figure'));
clc

%% Run Simulation
sim = drv.Simulation();

% Setup mesh result options:
sim.res.eid = true;
sim.res.nid = true;

% Setup nodal result options:
sim.res.scl = 1.0;
sim.res.dx  = false;  sim.res.dy  = false;  sim.res.dz  = false;
sim.res.rx  = false;  sim.res.ry  = false;  sim.res.rz  = false;

% Setup pption for smoothing element results:
sim.res.smooth = true;

% Setup element stress result options:
sim.res.sxx = true;   sim.res.syy = true;   sim.res.szz = false;
sim.res.txy = true;   sim.res.txz = false;  sim.res.tyz = false;
sim.res.s1  = true;   sim.res.s2  = true;   sim.res.s3  = false;
sim.res.taumax = true;

% Setup element internal forces result options:
sim.res.mxx = false;   sim.res.myy = false;
sim.res.qxz = false;   sim.res.qyz = false;
sim.res.m1  = false;   sim.res.m2  = false;
sim.res.tormax = false;

% Run simulation
sim.runAll();
