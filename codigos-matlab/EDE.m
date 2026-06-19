%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Metodo: Euler-Maruyama + Monte Carlo

close all; clear all; clc; tic;

%% Parametros del Modelo

% Parametros del Modelo (bajo P)
r = 0.05;
kappa = 2;
theta = 0.04;
sigma = 0.3;
rho = -0.5;
mu = r;
lambda = 0;

% Parametros de la Opcion
K = 80;
T = 1;
S0 = 100;
ups0 = 0.04;

% Parametros del Modelo (bajo Q)
kappa_tilde = kappa + lambda*sigma;
theta_tilde = (kappa*theta - rho*sigma*(mu - r))/kappa_tilde;

% Parametros de Monte Carlo
N_tray = 10000;
N_plot = 1000;
% N_rep = 1000;

% Numero de bloques para IC empirico 
B_bloques = 20;               
m = N_tray / B_bloques;

% Parametros de Discretizacion Temporal
N_t = 1000;
dt = T / N_t;
t = linspace(0, T, (N_t + 1));

%% Simulacion de Trayectorias

payoff_call = zeros(N_tray, 1);
payoff_put = zeros(N_tray, 1);
S = zeros((N_t + 1), 1);
ups = zeros((N_t + 1), 1);

S_plot = zeros((N_t + 1), N_plot);
ups_plot = zeros((N_t + 1), N_plot);

for n_tray = 1:N_tray
    S(1) = S0;
    ups(1) = ups0;
    dW1 = sqrt(dt)*randn(N_t, 1);
    dW2 = sqrt(dt)*randn(N_t, 1);

    for i = 1:N_t
        ups_iplus = max(ups(i), 0);

        S(i+1) = S(i) + r*S(i)*dt + sqrt(ups_iplus)*S(i)*dW1(i);
        ups(i+1) = ups_iplus + kappa_tilde*(theta_tilde - ups_iplus)*dt ...
            + sigma*rho*sqrt(ups_iplus)*dW1(i) ...
            + sigma*sqrt((1 - rho^2)*ups_iplus)*dW2(i);
    end

    payoff_call(n_tray) = max(S(end) - K, 0);
    payoff_put(n_tray) = max(K - S(end), 0);

    if n_tray <= N_plot
        S_plot(:, n_tray) = S;
        ups_plot(:, n_tray) = ups;
    end
end

descuento = exp(-r*T);
precio_call = descuento*mean(payoff_call);
precio_put = descuento*mean(payoff_put);

paridad_put_call = precio_call - precio_put - S0 + K*exp(-r*T);

tiempo = toc;

%% Intervalo de Confianza Numérico (teorico)

destip_call = std(payoff_call);
error_call = 1.96*destip_call/sqrt(N_tray);
IC_call = [precio_call - error_call, precio_call + error_call];

destip_put = std(payoff_put);
error_put = 1.96*destip_put/sqrt(N_tray);
IC_put = [precio_put - error_put, precio_put + error_put];

%% Intervalo de confianza empirico del estimador (por bloques, sin remuestreo)

precio_call_tray = descuento*payoff_call;
precio_put_tray  = descuento*payoff_put;

precio_call_bloques = zeros(B_bloques,1);
precio_put_bloques  = zeros(B_bloques,1);

for b = 1:B_bloques
    idx_ini = (b-1)*m + 1;
    idx_fin = b*m;
    precio_call_bloques(b) = mean(precio_call_tray(idx_ini:idx_fin));
    precio_put_bloques(b)  = mean(precio_put_tray(idx_ini:idx_fin));
end

IC_emp_call = prctile(precio_call_bloques,[2.5 97.5]);
IC_emp_put  = prctile(precio_put_bloques,[2.5 97.5]);
precio_call_mediana = median(precio_call_bloques);
precio_put_mediana  = median(precio_put_bloques);

%% Registro de Trayectorias de Muestra y Estadísticas

S_med = mean(S_plot, 2);
ups_med = mean(ups_plot, 2);
S_destip = std(S_plot, 0, 2);

S_T = zeros(N_tray, 1);

for n_tray = 1:N_tray
    if payoff_call(n_tray) > 0
        S_T(n_tray) = payoff_call(n_tray) + K;
    else
        S_T(n_tray) = K - payoff_put(n_tray);
    end
end

ups_T = ups_plot(end, :).';

%% Visualización de Figuras

figure(1); subplot(1,2,1); hold on; grid on;
for n_plot = 1:N_plot
    plot(t, S_plot(:, n_plot), 'Color', [[0.0 0.25 0.52], 1], 'LineWidth', 0.05);
    plot1 = plot(t, S_plot(:, n_plot), 'Color', [[0.0 0.25 0.52], 1], 'LineWidth', 0.05);
end
fill([t, fliplr(t)], [S_med.' + S_destip.', fliplr(S_med.' - S_destip.')], ...
    [0.6, 0.6, 0.6], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
fill1 = fill([t, fliplr(t)], [S_med.' + S_destip.', fliplr(S_med.' + S_destip.')], ...
    [0.6, 0.6, 0.6], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(t, S_med, 'Color', [0.1, 0.1, 0.1]);
plot2 = plot(t, S_med, 'Color', [0.1, 0.1, 0.1]);
yline(K, '--', 'Color', [0.2, 0.2, 0.2]);
yline1 = yline(K, '--', 'Color', [0.2, 0.2, 0.2]);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$S_t$', 'Interpreter', 'latex', 'FontSize', 18);
%title('Trayectorias $S_t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18); axis square;
legend([plot1, fill1, plot2, yline1], ...
    {'$S_t$','$[\bar{S}_t-\hat{\sigma}_{S_t},\ \bar{S}_t+\hat{\sigma}_{S_t}]$', '$\bar{S}_t$', '$K$'}, ...
    'Interpreter', 'latex', 'FontSize', 18);

figure(2); subplot(1,2,1); hold on; grid on;
for n_plot = 1:N_plot
    plot(t, ups_plot(:, n_plot), 'Color', [[0.0 0.25 0.52], 1], 'LineWidth', 0.05);
    plot1 = plot(t, ups_plot(:, n_plot), 'Color', [[0.0 0.25 0.52], 1], 'LineWidth', 0.05);
end
plot(t, ups_med, 'Color', [0.3 0.3 0.3]);
plot2 = plot(t, ups_med, 'Color', [0.3 0.3 0.3]);
yline(theta_tilde, '--', 'Color', [0.4 0.4 0.4]);
yline1 = yline(theta_tilde, '--', 'Color', [0.4 0.4 0.4]);
yline(ups0, ':', 'Color', [0.1 0.1 0.1]);
yline2 = yline(ups0, ':', 'Color', [0.1 0.1 0.1]);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\upsilon_t$', 'Interpreter', 'latex', 'FontSize', 18);
%title('Trayectorias $\\upsilon_t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);  axis square;
legend([plot1, plot2, yline1, yline2], ...
    {'$\upsilon_t$', '$\bar{\upsilon}_t$', '$\tilde{\theta}$', '$\upsilon_0$'}, ...
    'Interpreter', 'latex', 'FontSize', 18);

figure(1); subplot(1,2,2); hold on; grid on;
histogram(S_T, 80, 'Normalization', 'count', 'FaceColor', [0.1 0.35 0.62], 'EdgeColor', 'w');
xline(K, '--', 'Color', [0.3, 0.3, 0.3]);
xline(mean(S_T), '-', 'Color', [0.2, 0.2, 0.2]);
xline(median(S_T), ':', 'Color', [0.5, 0.5, 0.5]); axis square;
xlabel('$S_T$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\# S_T$', 'Interpreter', 'latex', 'FontSize', 18);
%title(sprintf('Distribucion $S_t$ - %d trayectorias', N_tray), ...
%'Interpreter', 'latex'); 
set(gca, 'FontSize', 18); axis square;
legend({'$\# S_T$', '$K$', '$\bar{S_T}$', '$\tilde{S_T}$'}, ...
    'Interpreter', 'latex', 'FontSize', 18);

figure(2); subplot(1,2,2); hold on; grid on;
a_cir = 2*kappa_tilde*theta_tilde/sigma^2;
b_cir = sigma^2/(2*kappa_tilde);
x_cir = linspace(0, max(ups_T)*1.5, 200);
pdf_cir = x_cir.^(a_cir-1).*exp(-x_cir/b_cir)/(b_cir^a_cir*gamma(a_cir));

histogram(ups_T, 30, 'Normalization', 'pdf', 'FaceColor', [0.1 0.35 0.62], 'EdgeColor', 'w');
plot(x_cir, pdf_cir, 'k', 'LineWidth', 2);
xline(theta_tilde, '--', 'Color', [0.3, 0.3, 0.3]); axis square;
xlabel('$\upsilon_T$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\# \upsilon_T$', 'Interpreter', 'latex', 'FontSize', 18);
%title(sprintf('Distribucion $\\upsilon_t$ - %d trayectorias', N_tray), ...
%'Interpreter', 'latex'); 
set(gca, 'FontSize', 18); axis square;
legend({'$\# \upsilon_T$', '$\Gamma(x; a_{\mathrm{CIR}}, b_{\mathrm{CIR}})$', ...
    '$\tilde{\theta}$'}, 'Interpreter', 'latex', 'FontSize', 18);

figure(3); subplot(1,2,1); hold on; grid on;
call_conv = descuento*cumsum(payoff_call)./(1:N_tray).';
semilogx(1:N_tray, call_conv, 'Color', [0.9, 0.4, 0.4]);
yline(precio_call, '--', 'Color', [0.4, 0.4, 0.4]);
xlabel('$n_{T}$', 'Interpreter', 'latex', 'FontSize', 18);
%ylabel('$\widehat{C}_{n_{T}}(0,S_0,\upsilon_0)$', 'Interpreter', 'latex', 'FontSize', 18);
%title('Convergencia del Estimador de Monte Carlo');
set(gca, 'FontSize', 18);  axis square;
legend({'$\widehat{C}^{\hspace{0.3mm}\textrm{EDE}}_{n_{T}}\hspace{0.1mm}(0,S_0,\upsilon_0)$ ', '$C_0$'}, ...
    'Interpreter', 'latex', 'FontSize', 18);

figure(3); subplot(1,2,2); hold on; grid on;
put_conv = descuento*cumsum(payoff_put)./(1:N_tray).';
semilogx(1:N_tray, put_conv, 'Color', [0.4, 0.8, 0.4]);
yline(precio_put, '--', 'Color', [0.4, 0.4, 0.4]);
xlabel('$n_{T}$', 'Interpreter', 'latex', 'FontSize', 18);
%ylabel('$\widehat{P}_{n_{T}}(0,S_0,\upsilon_0)$', 'Interpreter', 'latex', 'FontSize', 18);
%title('Convergencia del Estimador de Monte Carlo');
set(gca, 'FontSize', 18);  axis square;
legend({'$\widehat{P}^{\hspace{0.3mm}\textrm{EDE}}_{n_{T}}\hspace{0.1mm}(0,S_0,\upsilon_0)$ ', '$P_0$'}, ...
    'Interpreter', 'latex', 'FontSize', 18);

%% Mostrar Resultados

fprintf('==============================================================\n');
fprintf(' Valoracion de opciones europeas en el Modelo de Heston \n');
fprintf(' Metodo en EDEs: Euler-Maruyama + Monte Carlo \n');
fprintf('==============================================================\n');
fprintf('Precio Call: %.6f. (IC 95%%: [%.6f, %.6f])\n', precio_call, IC_call(1), IC_call(2));
fprintf('Precio Put: %.6f. (IC 95%%: [%.6f, %.6f])\n', precio_put, IC_put(1), IC_put(2));
%{
fprintf('Precio Call (empirico): media = %.6f, mediana = %.6f. (IC emp 95%%: [%.6f, %.6f])\n', ...
        mean(precio_call_bloques), precio_call_mediana, IC_emp_call(1), IC_emp_call(2));
fprintf('Precio Put  (empirico): media = %.6f, mediana = %.6f. (IC emp 95%%: [%.6f, %.6f])\n', ...
        mean(precio_put_bloques), precio_put_mediana, IC_emp_put(1), IC_emp_put(2));
%}
fprintf('Verificacion, Paridad Put-Call: %d.\n', round(abs(paridad_put_call)));
fprintf('Tiempo de Cómputo: %d.\n', tiempo);
fprintf('==============================================================\n');