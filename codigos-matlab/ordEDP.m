%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Metodo en EDP: Ordenes de Convergencia con Diferencias Finitas

close all; clear all; clc;

%% Parametros del Modelo

% Parametros del Modelo (bajo P)
r     = 0.05;
kappa = 2;
theta = 0.04;
sigma = 0.3;
rho   = -0.5;
mu    = r;
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

X_vec_ref = [linspace(50*exp(r*T), 75*exp(r*T), 20).'; ...
             linspace(85*exp(r*T), 160*exp(r*T), 30).'];
ups_vec_ref = linspace(c, 0.50, 20);

%% ORDEN ESPACIAL: OPCION 1 

n_sp1    = 5;
N_X0_1   = 40;
N_ups0_1 = 20;

V_sp1     = cell(1, n_sp1);
h_X_lv1   = zeros(1, n_sp1);
h_ups_lv1 = zeros(1, n_sp1);
N_tau_lv1 = zeros(1, n_sp1);

for lv = 1:n_sp1
    N_Xl   = N_X0_1*2^(lv - 1);
    N_upsl = N_ups0_1*2^(lv - 1);

    [H_lv, X_vec_lv, ups_vec_lv, N_tau_lv1(lv)] = ...
        solver_fd_full(N_Xl, N_upsl, b, c, d, T, sigma, ...
                       kappa_tilde, theta_tilde, rho, K);

    V_sp1{lv}     = interp_precio_2D(H_lv, X_vec_lv, ups_vec_lv, X_vec_ref, ups_vec_ref, r, T);
    h_X_lv1(lv)   = b/N_Xl;
    h_ups_lv1(lv) = (d - c)/N_upsl;
end

V_ref_sp1 = V_sp1{end};

n_err_sp1   = n_sp1 - 1;
err_sp1     = zeros(1, n_err_sp1);
h_X_err1    = zeros(1, n_err_sp1);
h_ups_err1  = zeros(1, n_err_sp1);
N_tau_err1  = zeros(1, n_err_sp1);

for k = 1:n_err_sp1
    d2           = (V_sp1{k} - V_ref_sp1).^2;
    err_sp1(k)   = sqrt(sum(d2(:))/numel(d2));
    h_X_err1(k)  = h_X_lv1(k);
    h_ups_err1(k)= h_ups_lv1(k);
    N_tau_err1(k)= N_tau_lv1(k);
end

ord_sp1 = log(err_sp1(1:end-1)./err_sp1(2:end)) ./ ...
          log(h_X_err1(1:end-1)./h_X_err1(2:end));

%% ORDEN ESPACIAL: OPCION 2 

n_sp2    = 5;
N_X0_2   = 64;
N_ups0_2 = 32;

N_X_list2   = zeros(1, n_sp2);
N_ups_list2 = zeros(1, n_sp2);

for lv = 1:n_sp2
    N_X_list2(lv)   = round(N_X0_2*(1.5)^(lv - 1));
    N_ups_list2(lv) = round(N_ups0_2*(1.5)^(lv - 1));
end

V_sp2     = cell(1, n_sp2);
h_X_lv2   = zeros(1, n_sp2);
h_ups_lv2 = zeros(1, n_sp2);
N_tau_lv2 = zeros(1, n_sp2);

for lv = 1:n_sp2
    N_Xl   = N_X_list2(lv);
    N_upsl = N_ups_list2(lv);

    [H_lv, X_vec_lv, ups_vec_lv, N_tau_lv2(lv)] = ...
        solver_fd_full(N_Xl, N_upsl, b, c, d, T, sigma, ...
                       kappa_tilde, theta_tilde, rho, K);

    V_sp2{lv}     = interp_precio_2D(H_lv, X_vec_lv, ups_vec_lv, X_vec_ref, ups_vec_ref, r, T);
    h_X_lv2(lv)   = b/N_Xl;
    h_ups_lv2(lv) = (d - c)/N_upsl;
end

V_ref_sp2 = V_sp2{end};

n_err_sp2   = n_sp2 - 1;
err_sp2     = zeros(1, n_err_sp2);
h_X_err2    = zeros(1, n_err_sp2);
h_ups_err2  = zeros(1, n_err_sp2);
N_tau_err2  = zeros(1, n_err_sp2);

for k = 1:n_err_sp2
    d2           = (V_sp2{k} - V_ref_sp2).^2;
    err_sp2(k)   = sqrt(sum(d2(:))/numel(d2));
    h_X_err2(k)  = h_X_lv2(k);
    h_ups_err2(k)= h_ups_lv2(k);
    N_tau_err2(k)= N_tau_lv2(k);
end

ord_sp2 = log(err_sp2(1:end-1)./err_sp2(2:end)) ./ ...
          log(h_X_err2(1:end-1)./h_X_err2(2:end));

%% ORDEN ESPACIAL: OPCION 3 

n_sp3    = 5;
N_X0_3   = 81;      
N_ups0_3 = 32;     

N_X_list3   = zeros(1, n_sp3);
N_ups_list3 = zeros(1, n_sp3);

for lv = 1:n_sp3
    N_X_list3(lv)   = round(N_X0_3*(4/3)^(lv - 1));
    N_ups_list3(lv) = round(N_ups0_3*(4/3)^(lv - 1));
end

V_sp3     = cell(1, n_sp3);
h_X_lv3   = zeros(1, n_sp3);
h_ups_lv3 = zeros(1, n_sp3);
N_tau_lv3 = zeros(1, n_sp3);

for lv = 1:n_sp3
    N_Xl   = N_X_list3(lv);
    N_upsl = N_ups_list3(lv);

    [H_lv, X_vec_lv, ups_vec_lv, N_tau_lv3(lv)] = ...
        solver_fd_full(N_Xl, N_upsl, b, c, d, T, sigma, ...
                       kappa_tilde, theta_tilde, rho, K);

    V_sp3{lv}     = interp_precio_2D(H_lv, X_vec_lv, ups_vec_lv, X_vec_ref, ups_vec_ref, r, T);
    h_X_lv3(lv)   = b/N_Xl;
    h_ups_lv3(lv) = (d - c)/N_upsl;
end

V_ref_sp3 = V_sp3{end};

n_err_sp3   = n_sp3 - 1;
err_sp3     = zeros(1, n_err_sp3);
h_X_err3    = zeros(1, n_err_sp3);
h_ups_err3  = zeros(1, n_err_sp3);
N_tau_err3  = zeros(1, n_err_sp3);

for k = 1:n_err_sp3
    d2            = (V_sp3{k} - V_ref_sp3).^2;
    err_sp3(k)    = sqrt(sum(d2(:))/numel(d2));
    h_X_err3(k)   = h_X_lv3(k);
    h_ups_err3(k) = h_ups_lv3(k);
    N_tau_err3(k) = N_tau_lv3(k);
end

ord_sp3 = log(err_sp3(1:end-1)./err_sp3(2:end)) ./ ...
          log(h_X_err3(1:end-1)./h_X_err3(2:end));

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

X_vec_t   = linspace(0, b, N_X_t + 1).';
ups_vec_t = linspace(c, d, N_ups_t + 1);

for lv = 1:n_t
    h_taul = h_tau_base/2^(lv - 1);
    N_taul = ceil(T/h_taul);
    h_taul = T/N_taul;

    H_lv = solver_fd_htau(N_X_t, N_ups_t, b, c, d, T, h_taul, ...
                          sigma, kappa_tilde, theta_tilde, rho, K);

    V_t{lv}        = interp_precio_2D(H_lv, X_vec_t, ups_vec_t, X_vec_ref, ups_vec_ref, r, T);
    h_tau_list(lv) = h_taul;
    N_tau_list(lv) = N_taul;
end

n_err_t = n_t - 1;
err_t   = zeros(1, n_err_t);

for k = 1:n_err_t
    V_ref_t = 2*V_t{k+1} - V_t{k};
    d2_t    = (V_t{k} - V_ref_t).^2;
    err_t(k) = sqrt(sum(d2_t(:))/numel(d2_t));
end

ord_t = log(err_t(1:end-1)./err_t(2:end)) ./ ...
        log(h_tau_list(1:n_err_t-1)./h_tau_list(2:n_err_t));

%% Precio de Verificacion

[~, i_S] = min(abs(X_vec_ref - S0));
[~, i_v] = min(abs(ups_vec_ref - ups0));
precio_call = V_t{end}(i_S, i_v);
precio_put  = precio_call - S0 + K*exp(-r*T);

%% Creacion de Figuras

h_X_dense1   = logspace(log10(min(h_X_err1)), log10(max(h_X_err1)), 100);
h_ups_dense1 = logspace(log10(min(h_ups_err1)), log10(max(h_ups_err1)), 100);
h_X_dense2   = logspace(log10(min(h_X_err2)), log10(max(h_X_err2)), 100);
h_ups_dense2 = logspace(log10(min(h_ups_err2)), log10(max(h_ups_err2)), 100);
h_X_dense3   = logspace(log10(min(h_X_err3)), log10(max(h_X_err3)), 100);
h_ups_dense3 = logspace(log10(min(h_ups_err3)), log10(max(h_ups_err3)), 100);
h_tau_dense  = logspace(log10(min(h_tau_list(1:n_err_t))), ...
                        log10(max(h_tau_list(1:n_err_t))), 100);

figure(1);
%set(gcf, 'Position', [100, 100, 1300, 420]);

subplot(1, 3, 1); hold on;
p1 = loglog(h_X_err1, err_sp1, 'o', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 4);
p2 = loglog(h_X_dense1, err_sp1(1)*(h_X_dense1/h_X_err1(1)).^1, ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
p3 = loglog(h_X_dense1, err_sp1(1)*(h_X_dense1/h_X_err1(1)).^2, ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
p4 = loglog(h_X_err2, err_sp2, 's', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 5);
p5 = loglog(h_X_err3, err_sp3, '+', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 6);
xlabel('$h_X$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$e_{esp}^{(\cdot)}$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 16);
legend([p1, p4, p5, p2, p3], ...
       {'$e_{esp}^{(1)}$', '$e_{esp}^{(2)}$', '$e_{esp}^{(3)}$', ...
        '$\mathcal{O}(h_X)$', '$\mathcal{O}(h_X^2)$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest');
grid on; axis square; axis tight; hold off;

subplot(1, 3, 2); hold on;
p1 = loglog(h_ups_err1, err_sp1, 'o', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 4);
p2 = loglog(h_ups_dense1, err_sp1(1)*(h_ups_dense1/h_ups_err1(1)).^1, ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
p3 = loglog(h_ups_dense1, err_sp1(1)*(h_ups_dense1/h_ups_err1(1)).^2, ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
p4 = loglog(h_ups_err2, err_sp2, 's', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 5);
p5 = loglog(h_ups_err3, err_sp3, '+', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 6);
xlabel('$h_\upsilon$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$e_{esp}^{(\cdot)}$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 16);
legend([p1, p4, p5, p2, p3], ...
       {'$e_{esp}^{(1)}$', '$e_{esp}^{(2)}$', '$e_{esp}^{(3)}$', ...
        '$\mathcal{O}(h_\upsilon)$', '$\mathcal{O}(h_\upsilon^2)$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest');
grid on; axis square; axis tight; hold off;

subplot(1, 3, 3); hold on;
p1 = loglog(h_tau_list(1:n_err_t), err_t, 'o', 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 4);
p2 = loglog(h_tau_dense, err_t(1)*(h_tau_dense/h_tau_list(1)).^1, ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
p3 = loglog(h_tau_dense, err_t(1)*(h_tau_dense/h_tau_list(1)).^2, ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
xlabel('$h_{\tau}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$e_{\tau}$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 16);
legend([p1, p2, p3], {'$e_{\tau}$', '$\mathcal{O}(h_{\tau})$', '$\mathcal{O}(h_{\tau}^2)$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest');
grid on; axis square; axis tight; hold off;

%% Mostrar Resultados

fprintf('================================================================\n');
fprintf(' Valoracion de opciones europeas en el Modelo de Heston \n');
fprintf(' Metodo en EDP: Ordenes de Convergencia con Diferencias Finitas \n');
fprintf('================================================================\n');
fprintf('ORDEN ESPACIAL (OPCION 1):\n\n');
fprintf('%-6s %-8s %-10s %-16s %-10s\n', ...
        'N_X', 'N_ups', 'N_tau', 'e_esp', 'ord_sp.');
fprintf('%s\n', repmat('-', 1, 55));
for k = 1:n_err_sp1
    N_Xk   = N_X0_1*2^(k-1);
    N_upsk = N_ups0_1*2^(k-1);
    if k > 1
        fprintf('%-6d %-8d %-10d %-16.4e %-10.3f\n', ...
                N_Xk, N_upsk, N_tau_err1(k), err_sp1(k), ord_sp1(k-1));
    else
        fprintf('%-6d %-8d %-10d %-16.4e --\n', ...
                N_Xk, N_upsk, N_tau_err1(k), err_sp1(k));
    end
end
fprintf('\nMedia del orden de convergencia espacial (opcion 1): %.3f.\n', mean(ord_sp1));
fprintf('----------------------------------------------------------------\n');
fprintf('ORDEN ESPACIAL (OPCION 2):\n\n');
fprintf('%-6s %-8s %-10s %-16s %-10s\n', ...
        'N_X', 'N_ups', 'N_tau', 'e_esp', 'ord_sp.');
fprintf('%s\n', repmat('-', 1, 55));
for k = 1:n_err_sp2
    N_Xk   = N_X_list2(k);
    N_upsk = N_ups_list2(k);
    if k > 1
        fprintf('%-6d %-8d %-10d %-16.4e %-10.3f\n', ...
                N_Xk, N_upsk, N_tau_err2(k), err_sp2(k), ord_sp2(k-1));
    else
        fprintf('%-6d %-8d %-10d %-16.4e --\n', ...
                N_Xk, N_upsk, N_tau_err2(k), err_sp2(k));
    end
end
fprintf('\nMedia del orden de convergencia espacial (opcion 2): %.3f.\n', mean(ord_sp2));
fprintf('----------------------------------------------------------------\n');
fprintf('ORDEN ESPACIAL (OPCION 3):\n\n');
fprintf('%-6s %-8s %-10s %-16s %-10s\n', ...
        'N_X', 'N_ups', 'N_tau', 'e_esp', 'ord_sp.');
fprintf('%s\n', repmat('-', 1, 55));
for k = 1:n_err_sp3
    N_Xk   = N_X_list3(k);
    N_upsk = N_ups_list3(k);
    if k > 1
        fprintf('%-6d %-8d %-10d %-16.4e %-10.3f\n', ...
                N_Xk, N_upsk, N_tau_err3(k), err_sp3(k), ord_sp3(k-1));
    else
        fprintf('%-6d %-8d %-10d %-16.4e --\n', ...
                N_Xk, N_upsk, N_tau_err3(k), err_sp3(k));
    end
end
fprintf('\nMedia del orden de convergencia espacial (opcion 3): %.3f.\n', mean(ord_sp3));
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

function [H, X_vec, ups_vec, N_tau, h_tau] = solver_fd_full(N_X, N_ups, b, c, d, T, ...
                                                      sigma, kappa_tilde, theta_tilde, rho, K)
h_X   = b/N_X;
h_ups = (d - c)/N_ups;
h_tau_CFL = 1/(d*b^2/h_X^2 + sigma^2*d/h_ups^2);
N_tau = ceil(T/(0.9*h_tau_CFL));
h_tau = T/N_tau;

X_vec   = linspace(0, b, N_X + 1).';
ups_vec = linspace(c, d, N_ups + 1);
H       = cond_ini(X_vec, ups_vec, K);

for i_tau = 1:N_tau
    H = paso_euler_edp(H, X_vec, ups_vec, h_tau, h_X, h_ups, ...
                        sigma, kappa_tilde, theta_tilde, rho);
end
end

function H = solver_fd_htau(N_X, N_ups, b, c, d, T, h_tau, ...
                            sigma, kappa_tilde, theta_tilde, rho, K)
N_tau = round(T/h_tau);
h_tau = T/N_tau;
h_X   = b/N_X;
h_ups = (d - c)/N_ups;

X_vec   = linspace(0, b, N_X + 1).';
ups_vec = linspace(c, d, N_ups + 1);
H       = cond_ini(X_vec, ups_vec, K);

for i_tau = 1:N_tau
    H = paso_euler_edp(H, X_vec, ups_vec, h_tau, h_X, h_ups, ...
                        sigma, kappa_tilde, theta_tilde, rho);
end
end

function H = cond_ini(X_vec, ups_vec, K)
N_ups = length(ups_vec) - 1;
H     = repmat(max(X_vec - K, 0), 1, N_ups + 1);
H(1, :)   = 0;
H(end, :) = max(X_vec(end) - K, 0);
end

function V_out = interp_precio_2D(H, X_vec, ups_vec, X_vec_ref, ups_vec_ref, r, T)
N_rX  = length(X_vec_ref);
N_rU  = length(ups_vec_ref);
N_X   = length(X_vec);
N_ups = length(ups_vec);
V_out = zeros(N_rX, N_rU);

for i_X = 1:N_rX
    Xi = X_vec_ref(i_X);
    i0 = find(X_vec <= Xi, 1, 'last');
    if isempty(i0), i0 = 1; end
    if i0 >= N_X, i0 = N_X - 1; end
    i1 = i0 + 1;
    tX = (Xi - X_vec(i0))/(X_vec(i1) - X_vec(i0));

    for i_ups = 1:N_rU
        upsi = ups_vec_ref(i_ups);
        j0 = find(ups_vec <= upsi, 1, 'last');
        if isempty(j0), j0 = 1; end
        if j0 >= N_ups, j0 = N_ups - 1; end
        j1 = j0 + 1;
        tu = (upsi - ups_vec(j0))/(ups_vec(j1) - ups_vec(j0));

        V_out(i_X, i_ups) = exp(-r*T)*( ...
            (1 - tX)*(1 - tu)*H(i0, j0) + tX*(1 - tu)*H(i1, j0) + ...
            (1 - tX)*tu*H(i0, j1)      + tX*tu*H(i1, j1));
    end
end
end

function H_new = paso_euler_edp(H, X_vec, ups_vec, h_tau, h_X, h_ups, ...
                                 sigma, kappa_tilde, theta_tilde, rho)
N_X   = length(X_vec) - 1;
N_ups = length(ups_vec) - 1;
H_new = H;

for i_X = 2:N_X
    XiX = X_vec(i_X);
    for i_ups = 2:N_ups
        upsiups = ups_vec(i_ups);

        a_iX_iups = (h_tau/(2*h_X^2))*upsiups*XiX^2;
        b_iX_iups = (h_tau*rho*sigma*upsiups*XiX)/(4*h_X*h_ups);
        aux1      = (h_tau*sigma^2*upsiups)/(2*h_ups^2);
        aux2      = (h_tau*kappa_tilde*(theta_tilde - upsiups))/(2*h_ups);
        c_iX_iups = aux1 - aux2;
        d_iX_iups = aux1 + aux2;
        e_iX_iups = 1 - (h_tau/h_X^2)*upsiups*XiX^2 ...
                      - (h_tau*sigma^2/h_ups^2)*upsiups;

        H_new(i_X, i_ups) = ...
            a_iX_iups*(H(i_X-1, i_ups) + H(i_X+1, i_ups)) + ...
            b_iX_iups*(H(i_X+1, i_ups+1) - H(i_X+1, i_ups-1) ...
                     - H(i_X-1, i_ups+1) + H(i_X-1, i_ups-1)) + ...
            c_iX_iups*H(i_X, i_ups-1) + ...
            d_iX_iups*H(i_X, i_ups+1) + ...
            e_iX_iups*H(i_X, i_ups);
    end
end

% Columna i_ups = 1
for i_X = 2:N_X
    XiX  = X_vec(i_X);
    upsc = ups_vec(1);

    a_iX_cups = (h_tau/(2*h_X^2))*upsc*XiX^2;
    aux3      = (h_tau*kappa_tilde*(theta_tilde - upsc))/h_ups;
    e_iX_cups = 1 - (h_tau/h_X^2)*upsc*XiX^2 - aux3;

    H_new(i_X, 1) = ...
        a_iX_cups*(H(i_X-1, 1) + H(i_X+1, 1)) + ...
        aux3*H(i_X, 2) + ...
        e_iX_cups*H(i_X, 1);
end

% Columna i_ups = N_ups + 1
for i_X = 2:N_X
    XiX  = X_vec(i_X);
    upsd = ups_vec(end);

    a_iX_dups = (h_tau/(2*h_X^2))*upsd*XiX^2;
    aux3      = (h_tau*kappa_tilde*(theta_tilde - upsd))/h_ups;
    e_iX_dups = 1 - (h_tau/h_X^2)*upsd*XiX^2 + aux3;

    H_new(i_X, end) = ...
        a_iX_dups*(H(i_X-1, end) + H(i_X+1, end)) - ...
        aux3*H(i_X, end-1) + ...
        e_iX_dups*H(i_X, end);
end
end