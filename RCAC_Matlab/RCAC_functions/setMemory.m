%% Memory allocation

x           = zeros(lx, steps);
u           = zeros(lu, steps);
y           = zeros(ly,steps);
z           = zeros(lz,steps);
z_hat       = zeros(lz,steps);

phi_filt    = zeros(lz,ltheta, steps);
phi         = zeros(lphi,steps);
uf          = zeros(lz,steps);
P           = zeros(ltheta, ltheta,steps);
theta       = zeros(ltheta,steps);
