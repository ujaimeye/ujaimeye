Lx = 5; Ly = 5; Lz = 5;
Nx = 100; Ny = 100; Nz = 100;
tf = .5;
eps = 0.05;
theta2 = 0.3;
Pd_a = 1;

f_func = @(x, y) exp(-2*(x.^2 + y.^2));
f_x = @(x,y) -4*x.*exp(-2*(x.^2 + y.^2));

% Solve
[T, x, y, z, runtime] = euler(Lx, Ly, Lz, Nx, Ny, Nz, tf, eps, theta2, Pd_a, f_func, 1);
[T2, x2, y2, z2, runtime2] = euler(Lx, Ly, Lz, Nx, Ny, Nz, tf, eps, theta2, Pd_a, f_func, .5);
[T3, x3, y3, z3, runtime3] = euler(2*Lx, 2*Ly, 2*Lz, Nx, Ny, Nz, tf, eps, theta2, Pd_a, f_func, 1);
%[T4, x4, y4, z4, runtime4] = euler(Lx, Ly, Lz, 2*Nx, 2*Ny, 2*Nz, tf, eps, theta2, Pd_a, f_func, 1);

x0 = 0.05; y0 = 0.05;
[~, ix] = min(abs(x - x0));
[~, iy] = min(abs(y - y0));
%Solutions for error analysis 
T = squeeze(T(ix, iy, :));
T2 = squeeze(T2(ix, iy, :));
T3 =squeeze(T3(ix, iy, :));
%T4 =squeeze(T4(ix, iy, :));

Edt = T-T2 ; 
%E_dis = T-T4;
Edomain = T-T3 ; 

figure(1)
clf reset
plot(z,Edt,'-r','LineWidth',2)
xlabel('z');ylabel('error(T-T2)');title('Error for cdt=1 and cdt=.5')

figure(2)
clf reset
plot(z,Edomain,'-r','LineWidth',2)
xlabel('z');ylabel('error(T-T3)');title('Error for increased domain')


