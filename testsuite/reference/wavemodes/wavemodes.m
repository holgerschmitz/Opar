Ey = h5read('Ey.h5','/data');
By = h5read('By.h5','/data');
ey = fftshift(fft2(Ey'));
by = fftshift(fft2(By'));
Lx = 10;
T = 100;
s=size(Ey);
sx = s(1);
st = s(2);

delta_k = pi/Lx;
delta_om = pi/T;
klo = (-0.5-sx/2)*delta_k;
khi = (sx/2-1.5) *delta_k;
omlo = (-st/2)*delta_om;
omhi = (st/2-1) *delta_om;
klin = linspace(klo, khi, sx);
omlin = linspace(omlo, omhi, st);
z_coord = 0*omlin + 100;

[k_grid, om_grid] = meshgrid(klin, omlin);

% R and L mode
% (ck/om)^2 = 1 - (om_p/om)^2/(1 ± om_c/om)
% In normalised units: c = 1
% om_c = eB/m = 1
% om_p = sqrt(n e^2/ eps_0 m) = 1
% 

L_mode_k2 = omlin.^2 - 1./(1+1./omlin);
R_mode_k2 = omlin.^2 - 1./(1-1./omlin);

L_positive = L_mode_k2>0;
R_positive = R_mode_k2>0;

L_mode_k = sqrt(L_mode_k2);
R_mode_k = sqrt(R_mode_k2);

hold off;

surf(k_grid, om_grid, abs(ey),'EdgeColor','none','LineStyle','none'); 
view(0,90);
hold on;
plot3(L_mode_k(L_positive),omlin(L_positive),z_coord(L_positive),'.', 'linewidth', 2);
plot3(R_mode_k(R_positive),omlin(R_positive),z_coord(R_positive),'.', 'linewidth', 2);

xlim([-10,10]);
ylim([-10,10]);