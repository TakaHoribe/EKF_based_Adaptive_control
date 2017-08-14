function [A, B] = continuous_linear_dynamics_func(z_est, param)
% system includes position x

m = 1 / z_est(3);
b = param(2);
iz = 1 / z_est(4);
l = param(4);


A = [0, 1, 0
     0, -b/m, 0
     0, 0, 0];
 
B = [0, 0
     1/2/m, 1/2/m
     l/iz, -l/iz];