function x_next = discritized_dynamics_func(x, u, dt, param)

v = x(1);
w = x(2);
m_inv = x(3);
iz_inv = x(4);

u_r = u(1);
u_l = u(2);

% m = param(1);
b = param(2);
% iz = param(3);
l = param(4);


fx = [v - b * v * m_inv * dt
      w
      m_inv
      iz_inv];
  
gxu = [m_inv *  (u_r + u_l) / 2 * dt
       iz_inv * l * (u_r - u_l) * dt
       0
       0];

x_next = fx + gxu;
    
end