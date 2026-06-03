%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Ordenes de Convergencia del Metodo Euler-Maruyama + Monte Carlo

close all; clear all; clc;

%% Parametros del Modelo

% Parametros del Modelo (bajo P)
r      = 0.05;
kappa  = 2;
theta  = 0.04;
sigma  = 0.3;
rho    = -0.5;
mu     = r;
lambda = 0;

% Parametros de la Opcion
K    = 80;
T    = 1;
S0   = 100;
ups0 = 0.04;

% Parametros del Modelo (bajo Q)
kappa_tilde = kappa + lambda*sigma;
theta_tilde = (kappa*theta - rho*sigma*(mu - r))/kappa_tilde;

% Fecha de Valoracion y Vencimiento
fecha_valo = datenum(2026, 1, 1);
fecha_venc = datemnth(fecha_valo, round(T*12));

%% Precio de Referencia Semi-Analitico

call_ref = optByHestonNI(r, S0, fecha_valo, fecha_venc, 'call', K, ...
    ups0, theta_tilde, kappa_tilde, sigma, rho, ...
    'DividendYield', 0);

put_ref = optByHestonNI(r, S0, fecha_valo, fecha_venc, 'put', K, ...
    ups0, theta_tilde, kappa_tilde, sigma, rho, ...
    'DividendYield', 0);

%% Parametros del Estudio de Ordenes

expo_f   = 3:10; 
n_vec_f  = 2.^expo_f;
dt_vec_f = T./n_vec_f;
Lf       = numel(n_vec_f);

expo_d   = 3:10;
n_vec_d  = 2.^expo_d;
dt_vec_d = T./n_vec_d;
Ld       = numel(n_vec_d);

n_ref  = 2^12; 
dt_ref = T/n_ref;

M_conv     = 1000000;
SNR_umbral = 2;

rng(42);

%% Calculo de Errores Fuerte y Debil

error_fuerte = zeros(Lf, 1);
error_debil  = zeros(Ld, 1);
std_debil    = zeros(Ld, 1);

descuento = exp(-r*T);

blk = 50000;
nb  = M_conv/blk;

for k = 1:Lf
    n_k        = n_vec_f(k);
    dt_k       = dt_vec_f(k);
    ratio      = n_ref/n_k;
    calc_debil = (k <= Ld);

    acc_fuerte = 0;
    if calc_debil
        acc_diff  = 0;
        acc_diff2 = 0;
    end

    for b = 1:nb
        dW1_fino = sqrt(dt_ref)*randn(n_ref, blk);
        dW2_fino = sqrt(dt_ref)*randn(n_ref, blk);

        S_f   = S0*ones(1, blk);
        ups_f = ups0*ones(1, blk);

        for i = 1:n_ref
            ups_p = max(ups_f, 0);
            S_f   = S_f   + r*S_f.*dt_ref   + sqrt(ups_p).*S_f.*dW1_fino(i,:);
            ups_f = ups_p + kappa_tilde*(theta_tilde - ups_p)*dt_ref ...
                  + sigma*rho*sqrt(ups_p).*dW1_fino(i,:) ...
                  + sigma*sqrt(max(1 - rho^2, 0)*ups_p).*dW2_fino(i,:);
        end

        S_T_fino = S_f;

        dW1_gros = squeeze(sum(reshape(dW1_fino, ratio, n_k, blk), 1));
        dW2_gros = squeeze(sum(reshape(dW2_fino, ratio, n_k, blk), 1));

        S_g   = S0*ones(1, blk);
        ups_g = ups0*ones(1, blk);

        for i = 1:n_k
            ups_p = max(ups_g, 0);
            S_g   = S_g   + r*S_g.*dt_k   + sqrt(ups_p).*S_g.*dW1_gros(i,:);
            ups_g = ups_p + kappa_tilde*(theta_tilde - ups_p)*dt_k ...
                  + sigma*rho*sqrt(ups_p).*dW1_gros(i,:) ...
                  + sigma*sqrt(max(1 - rho^2, 0)*ups_p).*dW2_gros(i,:);
        end

        S_T_gros = S_g;

        acc_fuerte = acc_fuerte + sum(abs(S_T_gros - S_T_fino));

        if calc_debil
            d         = descuento*(max(S_T_gros - K, 0) - max(S_T_fino - K, 0));
            acc_diff  = acc_diff  + sum(d);
            acc_diff2 = acc_diff2 + sum(d.^2);
        end
    end

    error_fuerte(k) = acc_fuerte/M_conv;

    if calc_debil
        mu_d           = acc_diff/M_conv;
        var_d          = acc_diff2/M_conv - mu_d^2;
        error_debil(k) = abs(mu_d);
        std_debil(k)   = sqrt(max(var_d, 0)/M_conv);
    end
end

%% Ordenes Locales

snr_vec = error_debil./std_debil;

ord_local_f = NaN(Lf, 1);
for k = 2:Lf
    ord_local_f(k) = log(error_fuerte(k)/error_fuerte(k-1)) / ...
                     log(dt_vec_f(k)/dt_vec_f(k-1));
end

ord_local_d  = NaN(Ld, 1);
for k = 2:Ld
    ord_local_d(k) = log(error_debil(k)/error_debil(k-1)) / ...
                     log(dt_vec_d(k)/dt_vec_d(k-1));
end

%% Ordenes Globales por Regresion Log-Log y Calculo de Error Estandar

x_f = log(dt_vec_f(:));
y_f = log(error_fuerte(:));
coef_f  = polyfit(x_f, y_f, 1);
gamma_s = coef_f(1);
yfit_f = polyval(coef_f, x_f);
%SSE_f  = sum((y_f - yfit_f).^2);
%SST_f  = sum((x_f - mean(x_f)).^2);
%std_s  = sqrt(SSE_f / (Lf - 2)) / sqrt(SST_f);

x_d = log(dt_vec_d(:));
y_d = log(error_debil(:));
coef_w  = polyfit(x_d, y_d, 1);
gamma_w = coef_w(1);
yfit_d = polyval(coef_w, x_d);
%SSE_d  = sum((y_d - yfit_d).^2);
%SST_d  = sum((x_d - mean(x_d)).^2);
%std_w  = sqrt(SSE_d / (Ld - 2)) / sqrt(SST_d);

%% Visualizacion de Figuras

figure(1);

subplot(1, 2, 1); hold on; grid on;
log_dt_f = log(dt_vec_f(:));
log_es   = log(error_fuerte(:));
plot(log_dt_f, log_es, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 1.5);
plot(log_dt_f, yfit_f, ':', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);
xlabel('$\log (\Delta t)$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\log (e_{s})$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
legend({'$\log (e_{s})$', '$\mathcal{O}(\Delta t^{0.5})$'},'Interpreter', 'latex', 'FontSize', 18)
%txt_f = sprintf('Slope = %.4f, Std: %.4f', gamma_s, std_s);
%text(ax1, 0.05, 0.95, txt_f, 'Units', 'normalized', ...
%    'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'normal');
%title(ax1, 'Orden conv. fuerte', 'FontSize', 12, 'FontWeight', 'bold');
axis square; axis tight; hold off;

subplot(1, 2, 2); hold on; grid on;
log_dt_d = log(dt_vec_d(:));
log_ew   = log(error_debil(:));
plot(log_dt_d, log_ew, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 1.5);
plot(log_dt_d, yfit_d, ':', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);
xlabel('$\log (\Delta t)$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\log (e_{w})$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18);
legend({'$\log (e_{w})$', '$\mathcal{O}(\Delta t^{1})$'},'Interpreter', 'latex', 'FontSize', 18)
%txt_d = sprintf('Slope = %.4f, Std: %.4f', gamma_w, std_w);
%text(ax2, 0.05, 0.95, txt_d, 'Units', 'normalized', ...
%    'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'normal');
% title(ax2, 'Orden conv. débil', 'FontSize', 12, 'FontWeight', 'bold');
axis square; axis tight; hold off;

%% Mostrar Resultados

fprintf('===============================================================================\n');
fprintf('        Valoracion de opciones europeas en el Modelo de Heston                   \n');
fprintf('        Ordenes de Convergencia: Euler-Maruyama + Monte Carlo                    \n');
fprintf('===============================================================================\n');
sep  = repmat('-', 1, 79);
sep1 = repmat('-', 1, 79);
fprintf('%-6s %-10s %-12s %-9s %-12s %-9s %-8s\n', ...
    'n_k', 'dt_k', 'e_s', 'ord_s', 'e_w', 'ord_w', 'SNR');
fprintf('%s\n', sep1);
for k = 1:Lf
    ef_str = sprintf('%.3e', error_fuerte(k));
    if isnan(ord_local_f(k))
        of_str = '  -  ';
    else
        of_str = sprintf('%.4f', ord_local_f(k));
    end
    if k <= Ld
        ed_str = sprintf('%.3e', error_debil(k));
        if k == 1
            od_str = '  -  ';
        else
            od_str = sprintf('%.4f', ord_local_d(k));
        end
        snr_str = sprintf('%.2f', snr_vec(k));
    else
        ed_str  = '     -     ';
        od_str  = '    -    ';
        snr_str = '   -  ';
    end
    fprintf('%-6d %-10.6f %-12s %-9s %-12s %-9s %-8s\n', ...
        n_vec_f(k), dt_vec_f(k), ef_str, of_str, ed_str, od_str, snr_str);
end
fprintf('\n%s\n', sep1);
fprintf('Media del orden de convergencia fuerte: ord_s = %.4f \n', gamma_s);
fprintf('Media del orden de convergencia debil:  ord_w = %.4f \n', gamma_w);
fprintf('===============================================================================\n');
