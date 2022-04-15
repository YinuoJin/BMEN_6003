%% Q2
vmin = -40;
vmax = 60;
V = linspace(vmin, vmax, vmax-vmin+1);

n_v = 0.3177 % retrieved from Q1

alpha_n = 0.01 .* (10 - V) ./ ( exp((10-V)./10) - 1);
beta_n = 0.125 .* exp(-V/80);

alpha_n_new = 0.02 .* (10 - V) ./ ( exp((10-V)./10) - 1);
beta_n_new = 0.125 .* exp(-V/80);


n_inf = alpha_n ./ (alpha_n + beta_n);
n_inf_new = alpha_n_new ./ (alpha_n_new + beta_n_new);

dn_dv = n_inf(2:length(n_inf)) - n_inf(1:length(n_inf)-1);
dn_dv_new = n_inf_new(2:length(n_inf_new)) - n_inf_new(1:length(n_inf_new)-1);

V_k = -12*1e-3;
V_l = 10.6*1e-3;
g_k = 30*1e-3; % mho / cm^2
g_l = 0.3*1e-3; % mho / cm^2
C_m = 1e-6; % F / cm^2



%(1)
r = 1 / (g_k*4*n_v^3*(0 - V_k)*dn_dv(V==0)) * 1e-3;
r_new = 1 / (g_k*4*n_v^3*(0 - V_k)*dn_dv_new(V==0)) * 1e-3;
R = 1 / (g_l + g_k*n_v^4);

L = r / (alpha_n(V==0) + beta_n(V==0)) * 1e-3;
L_new = r_new / (alpha_n_new(V==0) + beta_n_new(V==0)) * 1e-3;

w_0 = sqrt(1 / (L * C_m));
Q = w_0 / (r/L + (1/(R*C_m)));


w_0_new = sqrt(1 / (L_new * C_m));
Q_new = w_0_new / (r_new/L + (1/(R*C_m)));

f_0 = 1 / w_0;
f_0_new = 1 / w_0_new;

disp([w_0, Q]);
disp([w_0_new, Q_new]);

% plot
figure(3)
subplot(1,2,1)
plot(V, alpha_n_new);
hold on;
plot(V, beta_n_new);
plot(V, n_inf_new);
legend('alpha_n', 'beta_n', 'n_inf');
xlabel("V = V_m-V_r (mV)");
ylabel("units 1/ms for alpha and beta, unitless for n_{inf}");

subplot(1,2,2)
plot(V(1:length(V)-1), dn_dv_new);
legend("dn_{inf} / d_V");
xlabel("V = V_m-V_r (mV)");
ylabel("units (1/mV)");

f = linspace(0, 200, 201);
w = 2*pi*f;

z1 = 1i.*w.*C_m;
z2 = 1/R;
z3 = 1./(r_new+1i.*w.*L_new);

Z = 1 ./ (z1+z2+z3);
Z_mag = abs(Z);
Z_phase = angle(Z) .* (180 / pi);

figure(4)
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




