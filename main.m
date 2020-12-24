%% FEMOOLab - Finite Element Method Object-Oriented Laboratory
%
%% Description
%
% This is the main script file of FEMOOLab.
% Run this script to select the input file with model and analysis data for
% performing the simulation.
%
% Multiple files can be selected, by holding down the Shift or Ctrl key and
% clicking on file names, to run the simulations sequentially.
%
% For more information, check the <readme.html README> file.
%
%% Plotting Options
% Set result plotting options with flags (true or false):

% Mesh labels:
opt.eid    = false;  % Plot element numbers (NOT IMPLEMENTED)
opt.nid    = false;  % Plot node numbers (NOT IMPLEMENTED)

% Nodal result:
opt.scl    = 1.0;    % Scale factor for deformed mesh (NOT IMPLEMENTED)
opt.dx     = false;  % Plot contour of displacements in X direction (NOT IMPLEMENTED)
opt.dy     = false;  % Plot contour of displacements in Y direction (NOT IMPLEMENTED)
opt.dz     = false;  % Plot contour of displacements in Z directiont (NOT IMPLEMENTED)
opt.rx     = false;  % Plot contour of rotations about X axis (NOT IMPLEMENTED)
opt.ry     = false;  % Plot contour of rotations about Y axis (NOT IMPLEMENTED)
opt.rz     = false;  % Plot contour of rotations about Z axis (NOT IMPLEMENTED)

% Smoothing:
opt.smooth = true;   % Smooth element results at common nodes

% Element stress result:
opt.sxx    = true;   % Plot contour of normal stresses in X direction
opt.syy    = true;   % Plot contour of normal stresses in Y direction
opt.szz    = false;  % Plot contour of normal stresses in Z direction
opt.txy    = true;   % Plot contour of XY shear stresses
opt.txz    = false;  % Plot contour of XZ shear stresses
opt.tyz    = false;  % Plot contour of YZ shear stresses
opt.s1     = true;   % Plot contour of principal stresses 1
opt.s2     = true;   % Plot contour of principal stresses 2
opt.s3     = false;  % Plot contour of principal stresses 3
opt.taumax = true;   % Plot contour of maximum shear stresses

% Element internal forces result:
opt.mxx    = false;  % Plot contour of moment about X direction (NOT IMPLEMENTED)
opt.myy    = false;  % Plot contour of moment about Y direction (NOT IMPLEMENTED)
opt.mxy    = false;  % Plot contour of torsion moment (NOT IMPLEMENTED)
opt.qxz    = false;  % Plot contour of XZ shear force (NOT IMPLEMENTED)
opt.qyz    = false;  % Plot contour of YZ shear force (NOT IMPLEMENTED)
opt.m1     = false;  % Plot contour of principal moment 1 (NOT IMPLEMENTED)
opt.m2     = false;  % Plot contour of principal moment 2 (NOT IMPLEMENTED)
opt.tormax = false;  % Plot contour of maximum torsion (NOT IMPLEMENTED)

%% Run Analysis
close(findall(0,'Type','figure')); clearvars -except opt; clc; 
drv.Simulation(opt);
