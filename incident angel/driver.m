% Driver

Lx = 5; Ly = 5; Lz = 5;
Nx = 1000; Ny = 1000; Nz = 1000;
tf = 1;
eps = 0.05;
theta2 = 0.3;
Pd_a = 1;
cdt = 1;


% Source function f(x, y)
f_func = @(x, y) exp(-2*(x.^2 + y.^2));
f_x = @(x,y) -4*x.*exp(-2*(x.^2 + y.^2));

% Solve
[T, x, y, z, runtime] = euler(Lx, Ly, Lz, Nx, Ny, Nz, tf, eps, theta2, Pd_a, f_func, cdt);

[T_leading,T_twoterm] = asymptotic(x,y,z,tf,f_func,f_x,theta2,eps);
x0 = 0.05; y0 = 0.05;
[~, ix] = min(abs(x - x0));
[~, iy] = min(abs(y - y0));
T_num = squeeze(T(ix, iy, :));



figure(1)
clf reset 
hold on
plot(z,T_leading,'-k','LineWidth',2)
plot(z,T_twoterm,'--r','LineWidth',2)
plot(z,T_num)