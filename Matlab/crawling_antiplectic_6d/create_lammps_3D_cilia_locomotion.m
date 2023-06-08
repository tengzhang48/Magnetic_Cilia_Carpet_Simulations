clc;    % Clear the command window.
fprintf('Beginning to run %s.m ...\n', mfilename);
close all;  % Close all figures (except those of imtool.)
%clearvars;
workspace;  % Make sure the workspace panel is showing.
addpath('./functions')  % All the functions stores in the subfolder

%%%%%%%%%%%%%%%%%%%%%  Parameter Setting Starts  %%%%%%%%%%%%%%%%%%%%%%%%%

%%% Specify the output file name only, the suffix is not needed %%%
filename = 'loco_20x1';  

%%% Specify the scale for demensionless units %%%
%%%      [M]	  [L]	   [T]	    V	    F	      rho	      k_visc       E	        g
units=[1.0E-06	1.0E-03	1.0E-03	1.0E+00	1.000E-03	1.000E+03	1.000E-03	1.000E+03	1.000E+03];

%%% Specify the simulation domain: change lengthes in the physical unit
Lx_b = 160e-3/units(2); % normlized by length unit
Ly_b = 8e-3/units(2); % normlized by length unit
Lz_b = 16e-3/units(2); % normlized by length unit
b_s = [0 Lx_b; 0 Ly_b; 0 Lz_b]; 

%%% Material properties: change the correspongding values in the physical units %%%
% Solid 1
mu = 6.208e4/units(8);  % normlized by stress unit      % shear modulus
lambda = 3.042e6/units(8);           % Lame constant
rho = 3000/units(6);  % normlized by density unit      % Density
gravity = -9.8/units(9); % normlized by acceleration      % Gravity

% Solid 2
mu2 = mu/3;  % normlized by stress unit      % shear modulus
lambda2 = lambda/3;           % Lame constant
rho2 = 1070/units(6);  % normlized by density unit      % Density
   
                                                         
%%% Parameters of generateBox3d(Lx,Ly,Lz,N1,N2,N3): generate a box with the corner (x_min,y_min,z_min) at (0,0,0) %%% 
% Box mesh dimension: change lengthes in the physical unit
Lx = 15350e-6/units(2); % normlized by length unit
Ly = 2500e-6/units(2); % normlized by length unit
Lz = 200e-6/units(2); % normlized by length unit
% Distritization: specify the number of elements along x,y,z direction
N1 = 307;
N2 = 50;
N3 = 10;

%%%%%%%%%%%%%%%%%%% Create solids using class solid %%%%%%%%%%%%%%%%%%%
%%% Create mesh by generateBox3d(Lx,Ly,Lz,N1,N2,N3) 

S1=solid(mu,rho,lambda,'locomotion_20x1_node.csv','locomotion_20x1_element.csv',1200,mu2,rho2,lambda2,1);  
S1.move(6,2,0);  % move the mesh along a vector

%%% Put all the solids you want to output in an array
Solid_array=[S1];

%%%%%%%%%%%%%%%%%%%%%  Parameter Setting Ends  %%%%%%%%%%%%%%%%%%%%%%%%%%%
M_total=assemble_model_with_magnetic_force(Solid_array);

write_lattice(filename, M_total ,b_s);
write_magnetic_forces(M_total);
write_gravity(M_total,gravity);

fprintf('Done running %s.m\n', mfilename);