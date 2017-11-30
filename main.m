%% MSE 380 Simulation
% Harrison Handley
% Michael Graves
% Richard Arthurs
% Thomas Krammer

%% Setup
clear all; close all;
addpath('nonlinear', 'linear')
params = params(); % these are global parameters used across various models. Change in params.m
simulation_time = 15;

%% Nonlinear Model
y0 = [0,0,0,0,0,0,0,0];
[t_nonlin,y_nonlin] = ode45(@(t,y) simulate_generation_nonlin(t,y,params), [0 simulation_time], y0);
[p_H2, p_O2] = gaspressure_nonlin(y_nonlin,params); 
plot_generation_nonlin(t_nonlin,y_nonlin,p_H2,p_O2);

%% Linear Model
[ t_lin, y_lin, moles_H2, moles_O2, p_H2, p_O2 ] = simulate_generation_linear(simulation_time,params);

