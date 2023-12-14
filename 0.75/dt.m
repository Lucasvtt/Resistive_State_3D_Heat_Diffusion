clear all;

dx = 0.4;
dy = 0.4;
dz = 0.4;
gamma = 10;
u = 5.79;

delta = 2/(1/dx^2+1/dy^2+1/dz^2);
delta_t_1 = u*(delta/4)*(1/sqrt(1+gamma^2));

dx = 0.4;
dy = 0.4;
dz = 0.2;
gamma = 10;
u = 5.79;

delta = 2/(1/dx^2+1/dy^2+1/dz^2);
delta_t_2 = u*(delta/4)*(1/sqrt(1+gamma^2));

p = delta_t_2/delta_t_1