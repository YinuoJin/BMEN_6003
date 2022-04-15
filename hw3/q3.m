%% q3
% Declare consts
V_k = -12*1e-3;
V_na = 115*1e-3;
V_l = 10.6*1e-3;
g_k = 30*1e-3; % mho / cm^2
g_l = 0.3*1e-3; % mho / cm^2
g_na = 90*1e-3; % mho / cm^2
C_m = 1e-6; % F / cm^2
h_v = 0.79;

% Calculate variables
V = linspace(-5, 5, 11);

alpha_n = 0.01 .* (10 - V) ./ ( exp((10-V)./10) - 1);
beta_n = 0.125 .* exp(-V/80);

n_inf = alpha_n ./ (alpha_n + beta_n);
dn_dv = n_inf(2:length(n_inf)) - n_inf(1:length(n_inf)-1);

alpha_m = 0.1 .* (25 - V) ./ ( exp((25-V)./10) - 1);
beta_m = 4 .* exp(-V/18);

m_inf = alpha_m ./ (alpha_m + beta_m);
dm_dv = m_inf(2:length(m_inf)) - m_inf(1:length(m_inf)-1);

n_v = n_inf(V==0);
m_v = m_inf(V==0);

R = 1 / (g_l + g_k*n_v^4 + g_na * m_v^3 * h_v);

r_n = 1 / (g_k*4*n_v^3*(0 - V_k)*dn_dv(V==0)) * 1e-3;
r_m = 1 / (g_na*3*m_v^2*h_v*(0 - V_na)*dm_dv(V==0)) * 1e-3;

L_n = r_n / ( alpha_n(V==0) + beta_n(V==0) ) * 1e-3;
L_m = 0; % Approximation

% (1). Calculate R'
R_prime = (R * r_m) / (R + r_m);

disp(R_prime);

% (2). Calculate w_0 & Q
w_0 = sqrt(1 / (L_n*C_m));
r = r_n;
Q = w_0 / ((r_n/L_n) + (1/R_prime*C_m));
disp([w_0, Q]);

% (3). 
f = linspace(0, 200, 201);
w = 2*pi*f;

z1 = 1i.*w.*C_m;
z2 = 1/R_prime;
z3 = 1./(r+1i.*w.*L_n);

Z = 1 ./ (z1+z2+z3);
Z_mag = abs(Z);
Z_phase = angle(Z) .* (180 / pi);

figure(6)
subplot(2,1,1)
plot(f, Z_mag);
ylim([0, 4000]);
grid on;
xlabel("frequency (Hz)");
ylabel("Magnitude of Z (ohms)");

subplot(2,1,2)
plot(f, Z_phase);
grid on;
xlabel("frequency (Hz)");
ylabel("Phase of Z (degree)");





