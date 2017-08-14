function [F, H] = linearized_dynamics_func(x, u, dt, param)

v = x(1);
% w = x(2);
m_inv = x(3);
% iz_inv = x(4);

u_r = u(1);
u_l = u(2);

% m = param(1);
b = param(2);
% iz = param(3);
l = param(4);

F = [ 1 - b * m_inv * dt, 0, (-b * v + (u_r + u_l) / 2) * dt, 0;
      0, 1, 0, l * (u_r - u_l) * dt
      0, 0, 1, 0
      0, 0, 0, 1];
H = [1, 0, 0, 0
     0, 1, 0, 0];
    
end