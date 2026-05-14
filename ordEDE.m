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

expo_d   = 3:6;
n_vec_d  = 2.^expo_d;
dt_vec_d = T./n_vec_d;
Ld       = numel(n_vec_d);

n_ref  = 2^12;
dt_ref = T/n_ref;

M_conv     = 2000000;
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
flag_local_d = blanks(Ld);

for k = 1:Ld
    if snr_vec(k) >= SNR_umbral
        flag_local_d(k) = ' ';
    else
        flag_local_d(k) = '*';
    end
end

for k = 2:Ld
    ambos_limpios = (snr_vec(k) >= SNR_umbral) && (snr_vec(k-1) >= SNR_umbral);
    monotono      = (error_debil(k) < error_debil(k-1));
    
    if ambos_limpios && monotono
        ord_local_d(k) = log(error_debil(k)/error_debil(k-1)) / ...
                         log(dt_vec_d(k)/dt_vec_d(k-1));
    elseif ambos_limpios && ~monotono
        flag_local_d(k) = 'M';
    end
end

%% Ordenes Globales por Regresion Log-Log

coef_f  = polyfit(log(dt_vec_f(:)), log(error_fuerte), 1);
gamma_s = coef_f(1);

idx_reg_d      = false(Ld, 1);
primero_limpio = find(snr_vec >= SNR_umbral, 1);
if ~isempty(primero_limpio)
    idx_reg_d(primero_limpio) = true;
end
for k = 2:Ld
    if snr_vec(k) >= SNR_umbral && error_debil(k) < error_debil(k-1)
        idx_reg_d(k) = true;
    end
end

n_reg_d = sum(idx_reg_d);
if n_reg_d >= 2
    coef_w     = polyfit(log(dt_vec_d(idx_reg_d)'), log(error_debil(idx_reg_d)), 1);
    gamma_w    = coef_w(1);
    gamma_w_ok = true;
else
    coef_w     = [NaN NaN];
    gamma_w    = NaN;
    gamma_w_ok = false;
end

%% Visualizacion de Figuras

%dt_fit_f = exp(linspace(log(dt_vec_f(end))*0.75, log(dt_vec_f(1))*1.1, 200));
%dt_fit_d = exp(linspace(log(dt_vec_d(end))*0.75, log(dt_vec_d(1))*1.1, 200));
dt_fit_f = linspace(0.001, 0.125, 200);
dt_fit_d = linspace(0.03, 0.125, 200);

ref_half = error_fuerte(1) * (dt_fit_f / dt_vec_f(1)).^0.5;

idx_d_ok = idx_reg_d & (error_debil > 0);

col_fuerte = [0.1 0.1 0.1];
col_debil  = [0.1 0.1 0.1];

figure(1);
ax1 = subplot(1, 2, 1);
hold(ax1, 'on');
hp1 = loglog(ax1, dt_vec_f, error_fuerte, '-', ...
    'Color', col_fuerte, 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', col_fuerte);
hp2 = loglog(ax1, dt_fit_f, ref_half, 'k-.', 'LineWidth', 1.2);
set(ax1, 'XDir', 'reverse', 'XGrid', 'on', 'YGrid', 'on', ...
    'GridAlpha', 0.35, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
xlabel(ax1, '$\Delta t$',   'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax1, '$e_s$', 'Interpreter', 'latex', 'FontSize', 12);
legend(ax1, [hp1, hp2], ...
    {'$e_s$', '$\mathcal{O}(\Delta t^{0.5})$'}, ...
    'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northwest');
axis square; axis tight; hold(ax1, 'off');

ax2 = subplot(1, 2, 2);
hold(ax2, 'on');
leg_handles_d = [];
leg_labels_d  = {};
if any(idx_d_ok)
    hd1 = loglog(ax2, dt_vec_d(idx_d_ok), error_debil(idx_d_ok), '-', ...
        'Color', col_debil, 'LineWidth', 1.8, ...
        'MarkerSize', 8, 'MarkerFaceColor', col_debil);
    leg_handles_d(end+1) = hd1;
    leg_labels_d{end+1}  = '$e_w$';
end
if ~isempty(primero_limpio)
    ref_one = error_debil(primero_limpio) * ...
              (dt_fit_d / dt_vec_d(primero_limpio)).^1.0;
    hd2 = loglog(ax2, dt_fit_d, ref_one, 'k:', 'LineWidth', 1.2);
    leg_handles_d(end+1) = hd2;
    leg_labels_d{end+1}  = '$\mathcal{O}(\Delta t^{1})$';
end
set(ax2, 'XDir', 'reverse', 'XGrid', 'on', 'YGrid', 'on', ...
    'GridAlpha', 0.35, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
xlabel(ax2, '$\Delta t$',  'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax2, '$e_w$', 'Interpreter', 'latex', 'FontSize', 12);
if ~isempty(leg_handles_d)
    legend(ax2, leg_handles_d, leg_labels_d, ...
        'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northwest');
end
axis square; axis tight; hold(ax2, 'off');

%% Mostrar Resultados

fprintf('===============================================================================\n');
fprintf('        Valoracion de opciones europeas en el Modelo de Heston                   \n');
fprintf('        Ordenes de Convergencia: Euler-Maruyama + Monte Carlo                    \n');
fprintf('===============================================================================\n');
%fprintf('Precio de referencia semi-analitico (Call):  %.6f\n', call_ref);
%fprintf('Precio de referencia semi-analitico (Put):   %.6f\n', put_ref);
%fprintf('Trayectorias por nivel: %d\n', M_conv);

sep = repmat('-', 1, 79);
sep1 = repmat('-', 1, 70);
%fprintf('%s\n', sep);
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
        fl     = flag_local_d(k);
        if isnan(ord_local_d(k))
            if k == 1
                od_str = '  -  ';
            elseif fl == '*'
                od_str = ' (*) ';
            elseif fl == 'M'
                od_str = ' (**) ';
            else
                od_str = '  -  ';
            end
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

fprintf('%s\n', sep1);
fprintf(' (*) SNR < %g en algun nivel del par (orden debil local no fiable)\n', SNR_umbral);
fprintf(' (**) e_debil no decrece al refinar (el ruido MC domina sobre el sesgo)\n');
fprintf('%s\n', sep);
fprintf('Media del orden de convergencia fuerte: ord_s = %.4f\n', gamma_s);
if gamma_w_ok
    fprintf('Media del orden de convergencia debil:  ord_w = %.4f\n', gamma_w);
else
    fprintf('Orden global debil: no disponible\n');
end
fprintf('===============================================================================\n');