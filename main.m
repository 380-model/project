%% MSE 380 Simulation
% Harrison Handley
% Michael Graves
% Richard Arthurs
% Thomas Krammer

%% Setup
clear all; close all;
addpath('nonlinear', 'linear')
params = params(); % These are global parameters used across various models. Change in params.m
generation_time = 360; % (seconds) time to generate energy
consumption_time = 360; % (seconds) time to simulate energy consumption
dt = 0.01;
%% Nonlinear Model
y0 = [0,0,0,0,0,0,0,0];
[t_nonlin,y_nonlin] = ode45(@(t,y) simulate_generation_nonlin(t,y,params), [0 generation_time], y0); % Simulate energy generation + storage 
[p_H2_nonlin, p_O2_nonlin] = gaspressure_nonlin(y_nonlin,params); 
plot_generation_nonlin(t_nonlin,y_nonlin);

energyConsumptionODE(p_O2_nonlin, p_H2_nonlin, dt,consumption_time); % simulate energy consumption side

%% Linear Model
[ t_lin, y_lin, moles_H2_lin, moles_O2_lin, p_H2_lin, p_O2_lin ] = simulate_generation_linear(generation_time,params);
energyConsumption(p_O2_lin, p_H2_lin, dt, consumption_time);

