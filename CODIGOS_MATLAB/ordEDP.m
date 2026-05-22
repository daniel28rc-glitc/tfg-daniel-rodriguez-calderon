%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Metodo en EDP: Ordenes de Convergencia con Diferencias Finitas

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

%% Dominio Truncado

b = 240;
c = 0.01;
d = 1;

%% Malla de Referencia Comun

X_ref   = [linspace(50*exp(r*T),  75*exp(r*T), 20).'; ...
           linspace(85*exp(r*T), 160*exp(r*T), 30).'];
ups_ref = linspace(c, 0.50, 20);

%% Orden Espacial

n_sp   = 5;
N_X0   = 40;
N_ups0 = 20;

V_sp     = cell(1, n_sp);
h_X_lv   = zeros(1, n_sp);
h_ups_lv = zeros(1, n_sp);
N_tau_lv = zeros(1, n_sp);

for lv = 1:n_sp
    N_Xl   = N_X0*2^(lv - 1);
    N_upsl = N_ups0*2^(lv - 1);

    [H_lv, Xs_lv, upss_lv, N_tau_lv(lv)] = ...
        solver_fd_full(N_Xl, N_upsl, b, c, d, T, sigma, kappa_tilde, theta_tilde, rho, K);

    V_sp{lv}     = interp_bil(H_lv, Xs_lv, upss_lv, X_ref, ups_ref, r, T);
    h_X_lv(lv)   = b/N_Xl;
    h_ups_lv(lv) = (d - c)/N_upsl;
end

n_err_sp  = n_sp - 1;
err_sp    = zeros(1, n_err_sp);
h_X_err   = zeros(1, n_err_sp);
h_ups_err = zeros(1, n_err_sp);
N_tau_err = zeros(1, n_err_sp);

for k = 1:n_err_sp
    V_ref_x      = 2*V_sp{k+1} - V_sp{k};
    d2           = (V_sp{k} - V_ref_x).^2;
    err_sp(k)    = sqrt(sum(d2(:))/numel(d2));
    h_X_err(k)   = h_X_lv(k);
    h_ups_err(k) = h_ups_lv(k);
    N_tau_err(k) = N_tau_lv(k);
end

ord_sp = log(err_sp(1:end-1)./err_sp(2:end)) ./ ...
         log(h_X_err(1:end-1)./h_X_err(2:end));

%% Orden Temporal

N_X_t   = 200;
N_ups_t = 100;
h_X_t   = b/N_X_t;
h_ups_t = (d - c)/N_ups_t;

h_tau_CFL  = 1/(d*b^2/h_X_t^2 + sigma^2*d/h_ups_t^2);
h_tau_base = 0.45*h_tau_CFL;
n_t        = 5;

V_t        = cell(1, n_t);
h_tau_list = zeros(1, n_t);
N_tau_list = zeros(1, n_t);

Xs_t   = linspace(0, b, N_X_t   + 1).';
upss_t = linspace(c, d, N_ups_t + 1);

for lv = 1:n_t
    h_taul = h_tau_base/2^(lv - 1);
    N_taul = ceil(T/h_taul);
    h_taul = T/N_taul;

    H_lv = solver_fd_htau(N_X_t, N_ups_t, b, c, d, T, h_taul, ...
        sigma, kappa_tilde, theta_tilde, rho, K);

    V_t{lv}        = interp_bil(H_lv, Xs_t, upss_t, X_ref, ups_ref, r, T);
    h_tau_list(lv) = h_taul;
    N_tau_list(lv) = N_taul;
end

n_err_t = n_t - 1;
err_t   = zeros(1, n_err_t);

for k = 1:n_err_t
    V_ref_t  = 2*V_t{k+1} - V_t{k};
    d2_t     = (V_t{k} - V_ref_t).^2;
    err_t(k) = sqrt(sum(d2_t(:))/numel(d2_t));
end

ord_t = log(err_t(1:end-1)./err_t(2:end)) ./ ...
        log(h_tau_list(1:n_err_t-1)./h_tau_list(2:n_err_t));

%% Precio de Verificacion (malla espacial mas fina)

[~, i_S] = min(abs(X_ref - S0));
[~, i_v] = min(abs(ups_ref - ups0));
precio_call = V_t{end}(i_S, i_v);
precio_put  = precio_call - S0 + K*exp(-r*T);

%% Creacion de Figuras

h_X_dense   = logspace(log10(min(h_X_err)),   log10(max(h_X_err)),   100);
h_ups_dense = logspace(log10(min(h_ups_err)), log10(max(h_ups_err)), 100);
h_tau_dense = logspace(log10(min(h_tau_list(1:n_err_t))), log10(max(h_tau_list(1:n_err_t))), 100);

figure(1);
set(gcf, 'Position', [100, 100, 1300, 420]);

subplot(1, 3, 1); hold on;
p1 = loglog(h_X_err, err_sp, '-', 'Color', [0.1 0.1 0.1], ...
    'MarkerFaceColor', [0.1 0.1 0.1], 'LineWidth', 1.8, 'MarkerSize', 7);
p2 = loglog(h_X_dense, err_sp(1)*(h_X_dense/h_X_err(1)).^1, ...
    'k--', 'LineWidth', 1.2);
p3 = loglog(h_X_dense, err_sp(1)*(h_X_dense/h_X_err(1)).^2, ...
    '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
xlabel('$h_X$', 'Interpreter', 'latex');
ylabel('$e_{esp}$', 'Interpreter', 'latex');
%title('Convergencia Espacial ($h_X$)', 'Interpreter', 'latex');
legend([p1, p2, p3], {'$e_{esp}$', '$\mathcal{O}(h_X)$', '$\mathcal{O}(h_X^2)$'}, ...
    'Interpreter', 'latex', 'FontSize', 10, 'Location', 'southwest');
grid on; axis square; hold off;

subplot(1, 3, 2); hold on;
p1 = loglog(h_ups_err, err_sp, '-', 'Color', [0.1 0.1 0.1], ...
    'MarkerFaceColor', [0.1 0.1 0.1], 'LineWidth', 1.8, 'MarkerSize', 7);
p2 = loglog(h_ups_dense, err_sp(1)*(h_ups_dense/h_ups_err(1)).^1, ...
    'k--', 'LineWidth', 1.2);
p3 = loglog(h_ups_dense, err_sp(1)*(h_ups_dense/h_ups_err(1)).^2, ...
    '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
xlabel('$h_\upsilon$', 'Interpreter', 'latex');
ylabel('$e_{esp}$', 'Interpreter', 'latex');
%title('Convergencia Espacial ($h_\upsilon$)', 'Interpreter', 'latex');
legend([p1, p2, p3], {'$e_{esp}$', '$\mathcal{O}(h_\upsilon)$', '$\mathcal{O}(h_\upsilon^2)$'}, ...
    'Interpreter', 'latex', 'FontSize', 10, 'Location', 'southwest');
grid on; axis square; hold off;

subplot(1, 3, 3); hold on;
p1 = loglog(h_tau_list(1:n_err_t), err_t, '-', 'Color', [0.1 0.1 0.1], ...
    'MarkerFaceColor', [0.1 0.1 0.1], 'LineWidth', 1.8, 'MarkerSize', 7);
p2 = loglog(h_tau_dense, err_t(1)*(h_tau_dense/h_tau_list(1)).^1, ...
    'k--', 'LineWidth', 1.2);
p3 = loglog(h_tau_dense, err_t(1)*(h_tau_dense/h_tau_list(1)).^2, ...
    '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
xlabel('$\Delta\tau$', 'Interpreter', 'latex');
ylabel('$e_t$', 'Interpreter', 'latex');
%title('Convergencia Temporal', 'Interpreter', 'latex');
legend([p1, p2, p3], {'$e_t$', '$\mathcal{O}(\Delta\tau)$', '$\mathcal{O}(\Delta\tau^2)$'}, ...
    'Interpreter', 'latex', 'FontSize', 10, 'Location', 'southwest');
grid on; axis square; hold off;

%% Mostrar Resultados

fprintf('================================================================\n');
fprintf(' Valoracion de opciones europeas en el Modelo de Heston \n');
fprintf(' Metodo en EDP: Ordenes de Convergencia con Diferencias Finitas \n');
fprintf('================================================================\n');
%fprintf('Precio de la opcion (malla espacial mas fina):\n');
%fprintf('  Precio Call: %f.\n', precio_call);
%fprintf('  Precio Put:  %f.\n', precio_put);
%fprintf('----------------------------------------------------------------\n');
fprintf('ORDEN ESPACIAL:\n\n');
fprintf('%-6s %-8s %-10s %-16s %-10s\n', ...
    'N_X', 'N_ups', 'N_tau', 'e_esp', 'ord_sp.');
fprintf('%s\n', repmat('-', 1, 55));
for k = 1:n_err_sp
    if k > 1
        fprintf('%-6d %-8d %-10d %-16.4e %-10.3f\n', ...
            N_X0*2^(k-1), N_ups0*2^(k-1), N_tau_err(k), err_sp(k), ord_sp(k-1));
    else
        fprintf('%-6d %-8d %-10d %-16.4e --\n', ...
            N_X0*2^(k-1), N_ups0*2^(k-1), N_tau_err(k), err_sp(k));
    end
end
fprintf('\nMedia del orden de convergencia espacial: %.3f.\n', mean(ord_sp(2:end)));
fprintf('----------------------------------------------------------------\n');
fprintf('ORDEN TEMPORAL:\n\n');
fprintf('%-14s %-10s %-16s %-10s\n', ...
    'h_tau', 'N_tau', 'e_t', 'ord_t.');
fprintf('%s\n', repmat('-', 1, 55));
for k = 1:n_err_t
    if k > 1
        fprintf('%-14.4e %-10d %-16.4e %-10.3f\n', ...
            h_tau_list(k), N_tau_list(k), err_t(k), ord_t(k-1));
    else
        fprintf('%-14.4e %-10d %-16.4e --\n', ...
            h_tau_list(k), N_tau_list(k), err_t(k));
    end
end
fprintf('\nMedia del orden de convergencia temporal: %.3f.\n', mean(ord_t));
fprintf('================================================================\n');

%% Funciones Locales

function [H, Xs, upss, N_tau, h_tau] = solver_fd_full(N_X, N_ups, b, c, d, T, ...
    sigma, kappa_tilde, theta_tilde, rho, K)

h_X       = b/N_X;
h_ups     = (d - c)/N_ups;
h_tau_CFL = 1/(d*b^2/h_X^2 + sigma^2*d/h_ups^2);
N_tau     = ceil(T/(0.9*h_tau_CFL));
h_tau     = T/N_tau;

Xs   = linspace(0, b, N_X   + 1).';
upss = linspace(c, d, N_ups + 1);
H    = cond_ini(Xs, upss, K);

for i_tau = 1:N_tau
    H = euler_explicito(H, Xs, upss, h_tau, h_X, h_ups, sigma, kappa_tilde, theta_tilde, rho);
end
end

function H = solver_fd_htau(N_X, N_ups, b, c, d, T, h_tau, ...
    sigma, kappa_tilde, theta_tilde, rho, K)

N_tau = round(T/h_tau);
h_tau = T/N_tau;
h_X   = b/N_X;
h_ups = (d - c)/N_ups;
Xs    = linspace(0, b, N_X   + 1).';
upss  = linspace(c, d, N_ups + 1);
H     = cond_ini(Xs, upss, K);

for i_tau = 1:N_tau
    H = euler_explicito(H, Xs, upss, h_tau, h_X, h_ups, sigma, kappa_tilde, theta_tilde, rho);
end
end

function H = cond_ini(Xs, upss, K)
N_ups = length(upss) - 1;
H     = repmat(max(Xs - K, 0), 1, N_ups + 1);
H(1,   :) = 0;
H(end, :) = max(Xs(end) - K, 0);
end

function V_out = interp_bil(H, Xs, upss, X_ref, ups_ref, r, T)
N_rX  = length(X_ref);
N_rU  = length(ups_ref);
N_X   = length(Xs);
N_ups = length(upss);
V_out = zeros(N_rX, N_rU);

for i_X = 1:N_rX
    Xi = X_ref(i_X);
    i0 = find(Xs <= Xi, 1, 'last');
    if isempty(i0), i0 = 1; end
    if i0 >= N_X,   i0 = N_X - 1; end
    i1 = i0 + 1;
    tX = (Xi - Xs(i0))/(Xs(i1) - Xs(i0));

    for i_ups = 1:N_rU
        upsi = ups_ref(i_ups);
        j0 = find(upss <= upsi, 1, 'last');
        if isempty(j0), j0 = 1; end
        if j0 >= N_ups, j0 = N_ups - 1; end
        j1 = j0 + 1;
        tu = (upsi - upss(j0))/(upss(j1) - upss(j0));

        V_out(i_X, i_ups) = exp(-r*T)*( ...
            (1 - tX)*(1 - tu)*H(i0, j0) + tX*(1 - tu)*H(i1, j0) + ...
            (1 - tX)*      tu*H(i0, j1) + tX*       tu*H(i1, j1));
    end
end
end

function H_new = euler_explicito(H, Xs, upss, h_tau, h_X, h_ups, ...
    sigma, kappa_tilde, theta_tilde, rho)

N_X   = length(Xs)   - 1;
N_ups = length(upss) - 1;
H_new = H;

% Nodos Interiores
for i_X = 2:N_X
    XiX = Xs(i_X);
    
    for i_ups = 2:N_ups
        upsiups = upss(i_ups);

        a_iX_iups = (h_tau/(2*h_X^2))*upsiups*XiX^2;
        b_iX_iups = (h_tau*rho*sigma*upsiups*XiX)/(4*h_X*h_ups);
        aux1      = (h_tau*sigma^2*upsiups)/(2*h_ups^2);
        aux2      = (h_tau*kappa_tilde*(theta_tilde - upsiups))/(2*h_ups);
        c_iX_iups = aux1 - aux2;
        d_iX_iups = aux1 + aux2;
        e_iX_iups = 1 - (h_tau/h_X^2)*upsiups*XiX^2 ...
            - (h_tau*sigma^2/h_ups^2)*upsiups;

        H_new(i_X, i_ups) = a_iX_iups*(H(i_X-1, i_ups) + H(i_X+1, i_ups)) ...
            + b_iX_iups*(H(i_X+1, i_ups+1) - H(i_X+1, i_ups-1) ...
            - H(i_X-1, i_ups+1) + H(i_X-1, i_ups-1)) ...
            + c_iX_iups*H(i_X, i_ups-1) ...
            + d_iX_iups*H(i_X, i_ups+1) ...
            + e_iX_iups*H(i_X, i_ups);
    end
end

% Columna i_ups = 1 (ups = c)
for i_X = 2:N_X
    XiX  = Xs(i_X);
    upsc = upss(1);

    a_iX_cups = (h_tau/(2*h_X^2))*upsc*XiX^2;
    aux3      = (h_tau*kappa_tilde*(theta_tilde - upsc))/h_ups;
    e_iX_cups = 1 - (h_tau/h_X^2)*upsc*XiX^2 - aux3;

    H_new(i_X, 1) = a_iX_cups*(H(i_X-1, 1) + H(i_X+1, 1)) ...
        + aux3*H(i_X, 2) ...
        + e_iX_cups*H(i_X, 1);
end

% Columna i_ups = N_ups + 1 (ups = d)
for i_X = 2:N_X
    XiX  = Xs(i_X);
    upsd = upss(end);

    a_iX_dups = (h_tau/(2*h_X^2))*upsd*XiX^2;
    aux3      = (h_tau*kappa_tilde*(theta_tilde - upsd))/h_ups;
    e_iX_dups = 1 - (h_tau/h_X^2)*upsd*XiX^2 + aux3;

    H_new(i_X, end) = a_iX_dups*(H(i_X-1, end) + H(i_X+1, end)) ...
        - aux3*H(i_X, end-1) ...
        + e_iX_dups*H(i_X, end);
end
end
