% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% with stress free boundary conditions,
% zero initial conditions
% and a time dependent force source,
% in a structured undeformed grid.
%
% Version 1a:	domain = rectangular
%            	medium = general (heterogeneous)
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%
% JBR 4/2020: Add absorbing boundary conditions
%
addpath('./functions/')

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here wavefield output : ****
avd_out = 'a'; % a=accel, v=vel, d=disp

%**** Set here the parameters of the square box domain and mesh : ****
P = 16; %4; %8; % polynomial degree (spatial resolution)
%********

%**** Set here the parameters of the time solver : ****
% Tduration = 500; %90*3; %150; %90; % seconds

%**** Set here source properties : ****
Fstr = 'ricker';
Ff0 = 1/3; %0.3; %1/10; %1/30; %1/2; %1/10; %0.3; % fundamental frequency
