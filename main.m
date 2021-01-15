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
opt.eid    = false;  % Plot element numbers
opt.nid    = false;  % Plot node numbers
opt.gid    = false;  % Plot gauss points

% Mesh deformation:
opt.deform = true;   % Plot deformed mesh
opt.scl    = 0.0;    % Scale factor for deformed mesh (false or 0.0: scale automatically calculated)

% Nodal results:
opt.dx     = true;   % Plot contour of displacements in X direction
opt.dy     = true;   % Plot contour of displacements in Y direction
opt.dz     = true;   % Plot contour of displacements in Z directiont
opt.rx     = true;   % Plot contour of rotations about X axis
opt.ry     = true;   % Plot contour of rotations about Y axis
opt.rz     = true;   % Plot contour of rotations about Z axis
opt.temp   = true;   % Plot contour of temperature field

% Smoothing:
opt.smooth = true;   % Smooth element results at common nodes

% Element stress results:
opt.sxx    = true;   % Plot contour of normal stresses in X direction
opt.syy    = true;   % Plot contour of normal stresses in Y direction
opt.szz    = true;   % Plot contour of normal stresses in Z direction
opt.txy    = true;   % Plot contour of XY shear stresses
opt.txz    = true;   % Plot contour of XZ shear stresses
opt.tyz    = true;   % Plot contour of YZ shear stresses
opt.s1     = true;   % Plot contour of principal stresses 1
opt.s2     = true;   % Plot contour of principal stresses 2
opt.s3     = true;   % Plot contour of principal stresses 3
opt.taumax = true;   % Plot contour of maximum shear stresses

% Element internal forces results:
opt.qxz    = true;   % Plot contour of XZ shear force
opt.qyz    = true;   % Plot contour of YZ shear force 
opt.mxx    = true;   % Plot contour of moment about X direction
opt.myy    = true;   % Plot contour of moment about Y direction
opt.mxy    = true;   % Plot contour of torsion moment
opt.m1     = true;   % Plot contour of principal moment 1
opt.m2     = true;   % Plot contour of principal moment 2
opt.tormax = true;   % Plot contour of maximum torsion

% Element heat flux results:
opt.fxx    = true;   % Plot contour of heat fluxes in X direction
opt.fyy    = true;   % Plot contour of heat fluxes in Y direction
opt.fzz    = true;   % Plot contour of heat fluxes in Z direction
opt.fm     = true;   % Plot contour of heat flux module

% Element properties:
opt.pec    = true;   % Plot contour of Peclet numbers

% Other options:
opt.tol    = 1e-5;   % Tolerance for cleaning small result values and differences

%% Run Analysis
close(findall(0,'Type','figure')); clearvars -except opt; clc; 
drv.Simulation(opt);
