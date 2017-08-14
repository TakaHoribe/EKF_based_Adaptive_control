function dx = continuous_dynamics_func(x, u, param)

v = x(1);
% w = x(2);
% m_inv = x(3);
% iz_inv = x(4);

u_r = u(1);
u_l = u(2);

m = param(1);
b = param(2);
iz = param(3);
l = param(4);


fx = [-b * v / m
      0];

gxu = [1 / m / 2 * (u_r + u_l)
       1 / iz * l * (u_r - u_l)];

dx = fx + gxu;
    
end