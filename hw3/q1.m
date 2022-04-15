%% Q1

% plot 1
vmin = -40;
vmax = 60;
V = linspace(vmin, vmax, vmax-vmin+1);
alpha_n = 0.01 .* (10 - V) ./ ( exp((10-V)./10) - 1);
beta_n = 0.125 .* exp(-V/80);
n_inf = alpha_n ./ (alpha_n + beta_n);
dn_inf = n_inf(2:length(n_inf)) - n_inf(1:length(n_inf)-1);
dV = V(2:length(V)) - V(1:length(V)-1);

figure(1)
subplot(1,2,1)
plot(V, alpha_n);
hold on;
plot(V, beta_n);
plot(V, n_inf);
legend('alpha_n', 'beta_n', 'n_inf');
xlabel("V = V_m-V_r (mV)");
ylabel("units 1/ms for alpha and beta, unitless for n_{inf}");

subplot(1,2,2)
plot(V(1:length(V)-1), dn_inf ./ dV);
legend("dn_{inf} / d_V");
xlabel("V = V_m-V_r (mV)");
ylabel("units (1/mV)");


% plot 2
% declaring variables, convert to MKS units
V_k = -12*1e-3;
V_lxww = 10.6*1e-3;
g_k = 30*1e-3; % mho / cm^2
g_l = 0.3*1e-3; % mho / cm^2
C_m = 1e-6; % F / cm^2

% find variable condition at V = 0 
n_v = n_inf(V == 0);
dn_dv = dn_inf(V == 0);

alpha_n_v = alpha_n(V == 0);
beta_n_v = beta_n(V == 0);

r = 1 / (g_k*4*n_v^3*(0 - V_k)*dn_dv) * 1e-3;
L = r / (alpha_n_v + beta_n_v) * 1e-3;
R = 1 / (g_l + g_k*n_v^4);

f = linspace(0, 200, 201);
w = 2*pi*f;

z1 = 1i.*w.*C_m;
z2 = 1/R;
z3 = 1./(r+1i.*w.*L);

Z = 1 ./ (z1+z2+z3);
Z_mag = abs(Z);
Z_phase = angle(Z) .* (180 / pi);

figure(2)
subplot(2,1,1)
plot(f, Z_mag);
ylim([0, 1400]);
grid on;
xlabel("frequency (Hz)");
ylabel("Magnitude of Z (ohms)");

subplot(2,1,2)
plot(f, Z_phase);
grid on;
xlabel("frequency (Hz)");
ylabel("Phase of Z (degree)");

