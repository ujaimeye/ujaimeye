function  [T,x_vals,y_vals,z_vals,runtime] = euler(Lx, Ly, Lz, Nx, Ny, Nz, tf, eps, theta2, Pd_a, f_func, cdt)

if nargin < 12
    cdt = 1.0;
end

tic; 

x_vals = linspace(0, Lx, Nx)';
y_vals = linspace(0, Ly, Ny)';
z_vals = linspace(0, Lz, Nz)';

dx = x_vals(2) - x_vals(1);
dy = y_vals(2) - y_vals(1);
dz = z_vals(2) - z_vals(1);

[x, y, z] = ndgrid(x_vals, y_vals, z_vals);

% === Stability-based time step calculation ===
lambda = cos(theta2);
dt_max = 1 / ((2*eps)/dx^2 + (2*eps)/dy^2 + 2/dz^2);
Nt = ceil(tf / dt_max);
dt_max_adj = tf / Nt;
dt = cdt * dt_max_adj;

fprintf('Using dt = %.2e, Nt = %d, dt_max_adj = %.2e\n', dt, Nt, dt_max_adj);

% === Source term ===
x_shift = x + eps * z * tan(theta2);
F_source = Pd_a * f_func(x, y) .* (1/lambda) .* exp(-z / lambda);

% === Initial condition ===
T = zeros(Nx, Ny, Nz);

% === Time-stepping ===
for n = 1:Nt
    T_old = T;

    d2Tdx2 = (T_old(3:end, 2:end-1, 2:end-1) - 2*T_old(2:end-1, 2:end-1, 2:end-1) + T_old(1:end-2, 2:end-1, 2:end-1)) / dx^2;
    d2Tdy2 = (T_old(2:end-1, 3:end, 2:end-1) - 2*T_old(2:end-1, 2:end-1, 2:end-1) + T_old(2:end-1, 1:end-2, 2:end-1)) / dy^2;
    d2Tdz2 = (T_old(2:end-1, 2:end-1, 3:end) - 2*T_old(2:end-1, 2:end-1, 2:end-1) + T_old(2:end-1, 2:end-1, 1:end-2)) / dz^2;

    S = F_source(2:end-1, 2:end-1, 2:end-1);

    T(2:end-1, 2:end-1, 2:end-1) = T_old(2:end-1, 2:end-1, 2:end-1) + dt * ...
        (eps^2 * (d2Tdx2 + d2Tdy2) + d2Tdz2 + S);

    % Neumann BC at z = 0
    T(:,:,1) = T(:,:,2);
end


runtime = toc;



end 

