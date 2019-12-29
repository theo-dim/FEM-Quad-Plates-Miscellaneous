%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stress.m - 1/12/16                                       %
% author: Tehila Stone | Theo Dimitrasopoulos              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; % removes all variables from the workspace.

function [sigma_x,sigma_y,tau_xy,tau_xz,tau,yz]=stress(E,v,epsilon_x,epsilon_y,gamma_xy,gamma_xz,gamma_yz)

L = zeros(5,1);
L = (E/(1-v^2))*[1 v 0 0 0; v 1 0 0 0; 0 0 0.5(1-v) 0 0; 0 0 0 0.5(1-v) 0; 0 0 0 0 0.5(1-v)]*[epsilon_x; epsilon_y; gamma_xy; gamma_xz; gamma_yz];
