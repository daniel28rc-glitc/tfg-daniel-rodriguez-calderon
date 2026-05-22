%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Comparacion de Metodos: Norma L^p(Q) via Feynman-Kac

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

%% Parametros de Simulacion

N_lp      = 3000;
t_eval    = 0.25;
n_SDE     = 500;
p_vec     = [1, 2];
Ng_vis    = 55;
M_mc      = 8000;
n_mc      = 200;
dt_mc     = T/n_mc;
descuento = exp(-r*T);

%% Nube bajo Q

rng(42);
dt_nube  = t_eval/n_SDE;
S_nube   = S0*ones(N_lp, 1);
ups_nube = ups0*ones(N_lp, 1);

for paso = 1:n_SDE
    dW1   = sqrt(dt_nube)*randn(N_lp, 1);
    dW2   = sqrt(dt_nube)*randn(N_lp, 1);
    ups_p = max(ups_nube, 0);
    S_nube   = S_nube   + r*S_nube.*dt_nube + sqrt(ups_p).*S_nube.*dW1;
    ups_nube = ups_p + kappa_tilde*(theta_tilde - ups_p)*dt_nube ...
             + sigma*rho*sqrt(ups_p).*dW1 ...
             + sigma*sqrt((1 - rho^2)*ups_p).*dW2;
end

ups_nube = max(ups_nube, 0);

%% Malla Visual

S_lo   = quantile(S_nube,   0.01); S_hi   = quantile(S_nube,   0.99);
ups_lo = quantile(ups_nube, 0.01); ups_hi = quantile(ups_nube, 0.99);

S_vis   = linspace(S_lo,   S_hi,   Ng_vis).';
ups_vis = linspace(ups_lo, ups_hi, Ng_vis).';
[SS_vis, UU_vis] = meshgrid(S_vis, ups_vis);
S_visv   = SS_vis(:);
ups_visv = max(UU_vis(:), 0);
N_vis    = length(S_visv);

%% EDP: Diferencias Finitas Explicitas

S_max_all   = max(max(S_nube),   max(S_visv));
ups_max_all = max(max(ups_nube), max(ups_visv));
b_dom = max(S_max_all*exp(r*T)*1.5, 3*K*exp(r*T));
c_dom = 0.001;
d_dom = max(ups_max_all*1.5, 1.0);

N_X   = 180;
N_ups = 80;
h_X   = b_dom/N_X;
h_ups = (d_dom - c_dom)/N_ups;

ma_CFL = d_dom*b_dom^2/h_X^2;
ms_CFL = sigma^2*d_dom/h_ups^2;
h_tau  = 0.85/(ma_CFL + ms_CFL);
N_tau  = ceil(T/h_tau);
h_tau  = T/N_tau;

Xs   = linspace(0,     b_dom, N_X   + 1).';
upss = linspace(c_dom, d_dom, N_ups + 1);

H_call = repmat(max(Xs - K, 0), 1, N_ups + 1);
H_put  = repmat(max(K - Xs, 0), 1, N_ups + 1);
H_call(1, :) = 0;  H_call(end, :) = max(b_dom - K, 0);
H_put(1,  :) = K;  H_put(end,  :) = 0;

for i_tau = 1:N_tau
    H_call = euler_explicito(H_call, Xs, upss, h_tau, h_X, h_ups, ...
                             sigma, kappa_tilde, theta_tilde, rho);
    H_put  = euler_explicito(H_put,  Xs, upss, h_tau, h_X, h_ups, ...
                             sigma, kappa_tilde, theta_tilde, rho);
end

F_EDP_call      = interpolacion_edp(S_nube,  ups_nube,  Xs, upss, H_call, r, T, b_dom, c_dom, d_dom);
F_EDP_put       = interpolacion_edp(S_nube,  ups_nube,  Xs, upss, H_put,  r, T, b_dom, c_dom, d_dom);
F_EDP_call_base = interpolacion_edp(S0,      ups0,      Xs, upss, H_call, r, T, b_dom, c_dom, d_dom);
F_EDP_put_base  = interpolacion_edp(S0,      ups0,      Xs, upss, H_put,  r, T, b_dom, c_dom, d_dom);
F_EDP_call_vis  = interpolacion_edp(S_visv,  ups_visv,  Xs, upss, H_call, r, T, b_dom, c_dom, d_dom);
F_EDP_put_vis   = interpolacion_edp(S_visv,  ups_visv,  Xs, upss, H_put,  r, T, b_dom, c_dom, d_dom);

%% EDE: Euler-Maruyama + Monte Carlo

rng(1234);
[F_EDE_call, F_EDE_put] = ede_lote(S_nube, ups_nube, M_mc, n_mc, dt_mc, descuento, ...
                                    K, kappa_tilde, theta_tilde, sigma, rho, r);

rng(5678);
payoff_call0 = zeros(M_mc, 1);
payoff_put0  = zeros(M_mc, 1);

for m = 1:M_mc
    St   = S0;
    upst = ups0;

    for s = 1:n_mc
        dW1   = sqrt(dt_mc)*randn;
        dW2   = sqrt(dt_mc)*randn;
        ups_p = max(upst, 0);
        St   = St   + r*St*dt_mc   + sqrt(ups_p)*St*dW1;
        upst = ups_p + kappa_tilde*(theta_tilde - ups_p)*dt_mc ...
             + sigma*rho*sqrt(ups_p)*dW1 ...
             + sigma*sqrt((1 - rho^2)*ups_p)*dW2;
    end

    payoff_call0(m) = max(St - K, 0);
    payoff_put0(m)  = max(K - St, 0);
end

F_EDE_call_base = descuento*mean(payoff_call0);
F_EDE_put_base  = descuento*mean(payoff_put0);

rng(3333);
[F_EDE_call_vis, F_EDE_put_vis] = ede_lote(S_visv, ups_visv, 4000, n_mc, dt_mc, descuento, ...
                                            K, kappa_tilde, theta_tilde, sigma, rho, r);

%% SA: Semi-Analitico 

F_SA_call_base = optByHestonNI(r, S0, fecha_valo, fecha_venc, 'call', K, ...
                               ups0, theta_tilde, kappa_tilde, sigma, rho, 'DividendYield', 0);
F_SA_put_base  = optByHestonNI(r, S0, fecha_valo, fecha_venc, 'put',  K, ...
                               ups0, theta_tilde, kappa_tilde, sigma, rho, 'DividendYield', 0);

F_SA_call = zeros(N_lp, 1);
F_SA_put  = zeros(N_lp, 1);

for k = 1:N_lp
    upsk = max(ups_nube(k), 1e-6);
    Sk   = max(S_nube(k),   0.01);
    try
        F_SA_call(k) = optByHestonNI(r, Sk, fecha_valo, fecha_venc, 'call', K, ...
                                     upsk, theta_tilde, kappa_tilde, sigma, rho, 'DividendYield', 0);
        F_SA_put(k)  = optByHestonNI(r, Sk, fecha_valo, fecha_venc, 'put',  K, ...
                                     upsk, theta_tilde, kappa_tilde, sigma, rho, 'DividendYield', 0);
    catch
        F_SA_call(k) = NaN;
        F_SA_put(k)  = NaN;
    end
end

F_SA_call_vis = zeros(N_vis, 1);
F_SA_put_vis  = zeros(N_vis, 1);

for k = 1:N_vis
    upsk = max(ups_visv(k), 1e-6);
    Sk   = max(S_visv(k),   0.01);
    try
        F_SA_call_vis(k) = optByHestonNI(r, Sk, fecha_valo, fecha_venc, 'call', K, ...
                                         upsk, theta_tilde, kappa_tilde, sigma, rho, 'DividendYield', 0);
        F_SA_put_vis(k)  = optByHestonNI(r, Sk, fecha_valo, fecha_venc, 'put',  K, ...
                                         upsk, theta_tilde, kappa_tilde, sigma, rho, 'DividendYield', 0);
    catch
        F_SA_call_vis(k) = NaN;
        F_SA_put_vis(k)  = NaN;
    end
end

%% Norma L^p(Q)

validos   = ~isnan(F_SA_call) & ~isnan(F_SA_put) & ...
            ~isnan(F_EDE_call) & ~isnan(F_EDP_call) & F_SA_call > 1e-8;
N_validos = sum(validos);

diff_EDE_SA_call  = abs(F_EDE_call(validos) - F_SA_call(validos));
diff_EDP_SA_call  = abs(F_EDP_call(validos) - F_SA_call(validos));
diff_EDE_EDP_call = abs(F_EDE_call(validos) - F_EDP_call(validos));
diff_EDE_SA_put   = abs(F_EDE_put(validos)  - F_SA_put(validos));
diff_EDP_SA_put   = abs(F_EDP_put(validos)  - F_SA_put(validos));
diff_EDE_EDP_put  = abs(F_EDE_put(validos)  - F_EDP_put(validos));

pares_call = {diff_EDE_SA_call,  'EDE vs SA'; ...
              diff_EDP_SA_call,  'EDP vs SA'; ...
              diff_EDE_EDP_call, 'EDE vs EDP'};
pares_put  = {diff_EDE_SA_put,   'EDE vs SA'; ...
              diff_EDP_SA_put,   'EDP vs SA'; ...
              diff_EDE_EDP_put,  'EDE vs EDP'};

Lp_call = zeros(3, length(p_vec));
Lp_put  = zeros(3, length(p_vec));

for i = 1:3
    for jp = 1:length(p_vec)
        Lp_call(i, jp) = mean(pares_call{i,1}.^p_vec(jp))^(1/p_vec(jp));
        Lp_put(i,  jp) = mean(pares_put{i,1} .^p_vec(jp))^(1/p_vec(jp));
    end
end

%% Curva L^p vs p

p_range       = 1:0.25:4;
Lp_curva_call = zeros(3, length(p_range));
Lp_curva_put  = zeros(3, length(p_range));
diffs_call    = {diff_EDE_SA_call, diff_EDP_SA_call, diff_EDE_EDP_call};
diffs_put     = {diff_EDE_SA_put,  diff_EDP_SA_put,  diff_EDE_EDP_put};
etiquetas     = {'EDE vs SA', 'EDP vs SA', 'EDE vs EDP'};

for i = 1:3
    for jp = 1:length(p_range)
        Lp_curva_call(i, jp) = mean(diffs_call{i}.^p_range(jp))^(1/p_range(jp));
        Lp_curva_put(i,  jp) = mean(diffs_put{i} .^p_range(jp))^(1/p_range(jp));
    end
end

colores = [0.1  0.1 0.1; ...
           0.1 0.1 0.1; ...
           0.1 0.1 0.1];

SA_malla_call  = reshape(F_SA_call_vis,  Ng_vis, Ng_vis);
EDE_malla_call = reshape(F_EDE_call_vis, Ng_vis, Ng_vis);
EDP_malla_call = reshape(F_EDP_call_vis, Ng_vis, Ng_vis);

SA_malla_put   = reshape(F_SA_put_vis,   Ng_vis, Ng_vis);
EDE_malla_put  = reshape(F_EDE_put_vis,  Ng_vis, Ng_vis);
EDP_malla_put  = reshape(F_EDP_put_vis,  Ng_vis, Ng_vis);

err_EDE_call_vis     = abs(EDE_malla_call - SA_malla_call);  err_EDE_call_vis(isnan(err_EDE_call_vis))     = 0;
err_EDP_call_vis     = abs(EDP_malla_call - SA_malla_call);  err_EDP_call_vis(isnan(err_EDP_call_vis))     = 0;
err_EDE_EDP_call_vis = abs(EDE_malla_call - EDP_malla_call); err_EDE_EDP_call_vis(isnan(err_EDE_EDP_call_vis)) = 0;
err_EDE_put_vis      = abs(EDE_malla_put  - SA_malla_put);   err_EDE_put_vis(isnan(err_EDE_put_vis))       = 0;
err_EDP_put_vis      = abs(EDP_malla_put  - SA_malla_put);   err_EDP_put_vis(isnan(err_EDP_put_vis))       = 0;
err_EDE_EDP_put_vis  = abs(EDE_malla_put  - EDP_malla_put);  err_EDE_EDP_put_vis(isnan(err_EDE_EDP_put_vis))  = 0;

[~, j_ups0] = min(abs(ups_vis - ups0));

figure(1);
subplot(1, 2, 1); hold on;
p1 = plot(S_vis, SA_malla_call(j_ups0,:).', '-',  'Color', [0.9, 0.4, 0.4], 'LineWidth', 2.0);
p2 = plot(S_vis, EDE_malla_call(j_ups0,:).','--', 'Color', [0.9, 0.4, 0.4], 'LineWidth', 1.8);
p3 = plot(S_vis, EDP_malla_call(j_ups0,:).',':', 'Color', [0.9, 0.4, 0.4], 'LineWidth', 1.8);
x1 = xline(K,  'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
x2 = xline(S0, 'k:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel('$S_t$', 'Interpreter', 'latex');
ylabel('$C_t$', 'Interpreter', 'latex');
%title(sprintf('Perfil Call ($\\upsilon_0 = %.4f$)', ups0), 'Interpreter', 'latex');
legend([p1, p2, p3, x1, x2], {'SA', 'EDE', 'EDP', '$K$', '$S_0$'}, 'Interpreter', 'latex', 'FontSize', 10);
grid on; axis square; hold off;

subplot(1, 2, 2); hold on;
p1 = plot(S_vis, SA_malla_put(j_ups0,:).', '-',  'Color', [0.4, 0.8, 0.4], 'LineWidth', 2.0);
p2 = plot(S_vis, EDE_malla_put(j_ups0,:).','--', 'Color', [0.4, 0.8, 0.4], 'LineWidth', 1.8);
p3 = plot(S_vis, EDP_malla_put(j_ups0,:).',':', 'Color', [0.4, 0.8, 0.4], 'LineWidth', 1.8);
x1 = xline(K,  'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
x2 = xline(S0, 'k:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel('$S_t$', 'Interpreter', 'latex');
ylabel('$P_t$', 'Interpreter', 'latex');
%title(sprintf('Perfil Put ($\\upsilon_0 = %.4f$)', ups0), 'Interpreter', 'latex');
legend([p1, p2, p3, x1, x2], {'SA', 'EDE', 'EDP', '$K$', '$S_0$'}, 'Interpreter', 'latex', 'FontSize', 10);
grid on; axis square; hold off;

figure(2); hold on;
p2 = scatter(S_nube, ups_nube, 6, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
p1 = plot(S0, ups0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('$S_t$', 'Interpreter', 'latex');
ylabel('$\upsilon_t$', 'Interpreter', 'latex');
%title('Nube Simulada bajo Medida Neutra', 'Interpreter', 'latex');
legend([p1, p2], {'$(S_0,\,\upsilon_0)$ ', 'Nube'}, 'Interpreter', 'latex', 'FontSize', 10);
grid on; axis square; hold off;

err_perfil_call = {err_EDE_call_vis(j_ups0,:).', ...
                   err_EDP_call_vis(j_ups0,:).', ...
                   abs(EDE_malla_call(j_ups0,:).' - EDP_malla_call(j_ups0,:).')};
err_malla_call  = {err_EDE_call_vis, err_EDP_call_vis, err_EDE_EDP_call_vis};
titulos_perfil  = {'$|EDE - SA|$', '$|EDP - SA|$', '$|EDE - EDP|$'};
titulos_mapa    = {'Mapa Error $|EDE - SA|$', 'Mapa Error $|EDP - SA|$', 'Mapa Error $|EDE - EDP|$'};
colores_cols    = [colores(1,:); colores(2,:); colores(3,:)];

figure(3);
for col = 1:3
    subplot(2, 3, col); hold on;
    plot(S_vis, err_perfil_call{col}, '-', 'Color', colores_cols(col,:), 'LineWidth', 2.0);
    xlabel('$S_t$', 'Interpreter', 'latex');
    %ylabel(titulos_perfil{col}, 'Interpreter', 'latex');
    %title(etiquetas{col}, 'Interpreter', 'latex');
    grid on; axis square; hold off; 

    subplot(2, 3, col + 3);
    clim_val = quantile(err_malla_call{col}(:), 0.97);
    imagesc(S_vis, ups_vis, err_malla_call{col}, [0, clim_val]);
    set(gca, 'YDir', 'normal'); colormap("gray"); colorbar;
    xlabel('$S_t$', 'Interpreter', 'latex');
    ylabel('$\upsilon_t$', 'Interpreter', 'latex');
    %title(titulos_mapa{col}, 'Interpreter', 'latex');
    hold on; axis square;
    xline(S0,   'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(ups0, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(S0, ups0, 'wo', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    hold off;
end

estilos = {'-', '-', '-'};

figure(4);
for col = 1:3
    subplot(1, 3, col); hold on;
    plot(p_range, Lp_curva_call(col,:), estilos{col}, ...
         'Color', colores_cols(col,:), 'LineWidth', 2, 'MarkerSize', 5, ...
         'DisplayName', 'Call');
    plot(p_range, Lp_curva_put(col,:), estilos{col}, ...
         'Color', colores_cols(col,:), 'LineWidth', 2, 'MarkerSize', 5, ...
         'LineStyle', '--', 'DisplayName', 'Put');
    xlabel('$p$', 'Interpreter', 'latex');
    %if col == 1
    %    ylabel('$\|EDE - SA\|_{L^p(Q)}$', 'Interpreter', 'latex');
    %elseif col == 2
    %    ylabel('$\|EDP - SA\|_{L^p(Q)}$', 'Interpreter', 'latex');
    %elseif col == 3
    %    ylabel('$\|EDE - EDP\|_{L^p(Q)}$', 'Interpreter', 'latex');
    %end
    %title(etiquetas{col}, 'Interpreter', 'latex');
    legend('Interpreter', 'latex', 'FontSize', 10);
    %xticks(1:0.5:4); 
    grid on; axis square; hold off;
end

%% Mostrar Resultados

fprintf('==============================================================\n');
fprintf('  Valoracion de opciones europeas en el Modelo de Heston     \n');
fprintf('  Comparacion de Metodos: EDE, EDP y SA - Norma L^p(Q)       \n');
fprintf('==============================================================\n');
%fprintf('Puntos validos de la nube: %d / %d.\n', N_validos, N_lp);
%fprintf('--------------------------------------------------------------\n');
fprintf('Precio puntual en (S0, ups0):\n\n');
fprintf('%-12s %12s %12s\n', 'Metodo', 'Call', 'Put');
fprintf('%s\n', repmat('-', 1, 40));
fprintf('%-12s %12.6f %12.6f\n', 'SA', F_SA_call_base, F_SA_put_base);
fprintf('%-12s %12.6f %12.6f\n', 'EDE',      F_EDE_call_base, F_EDE_put_base);
fprintf('%-12s %12.6f %12.6f\n', 'EDP',      F_EDP_call_base, F_EDP_put_base);
fprintf('--------------------------------------------------------------\n');
fprintf('CALL - Norma L^p(Q):\n\n');
fprintf('%-20s %12s %12s\n', 'Par de metodos', 'L^1(Q)', 'L^2(Q)');
fprintf('%s\n', repmat('-', 1, 48));
etiq = {'EDE vs SA', 'EDP vs SA', 'EDE vs EDP'};
for i = 1:3
    fprintf('%-20s %12.6f %12.6f\n', etiq{i}, Lp_call(i, 1), Lp_call(i, 2));
end
fprintf('--------------------------------------------------------------\n');
fprintf('PUT - Norma L^p(Q):\n\n');
fprintf('%-20s %12s %12s\n', 'Par de metodos', 'L^1(Q)', 'L^2(Q)');
fprintf('%s\n', repmat('-', 1, 48));
for i = 1:3
    fprintf('%-20s %12.6f %12.6f\n', etiq{i}, Lp_put(i, 1), Lp_put(i, 2));
end
fprintf('==============================================================\n');

%% Funciones Locales

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
        a = (h_tau/(2*h_X^2))*upsiups*XiX^2;
        b = (h_tau*rho*sigma*upsiups*XiX)/(4*h_X*h_ups);
        ts = (h_tau*sigma^2*upsiups)/(2*h_ups^2);
        td = (h_tau*kappa_tilde*(theta_tilde - upsiups))/(2*h_ups);
        c  = ts - td;
        d  = ts + td;
        e  = 1 - (h_tau/h_X^2)*upsiups*XiX^2 - (h_tau*sigma^2/h_ups^2)*upsiups;
        H_new(i_X, i_ups) = a*(H(i_X-1, i_ups) + H(i_X+1, i_ups)) ...
            + b*(H(i_X+1, i_ups+1) - H(i_X+1, i_ups-1) ...
               - H(i_X-1, i_ups+1) + H(i_X-1, i_ups-1)) ...
            + c*H(i_X, i_ups-1) + d*H(i_X, i_ups+1) + e*H(i_X, i_ups);
    end
end

% Columna i_ups = 1 (ups = c)
for i_X = 2:N_X
    XiX  = Xs(i_X); upsc = upss(1);
    a  = (h_tau/(2*h_X^2))*upsc*XiX^2;
    td = (h_tau*kappa_tilde*(theta_tilde - upsc))/h_ups;
    e  = 1 - (h_tau/h_X^2)*upsc*XiX^2 - td;
    H_new(i_X, 1) = a*(H(i_X-1, 1) + H(i_X+1, 1)) + td*H(i_X, 2) + e*H(i_X, 1);
end

% Columna i_ups = N_ups + 1 (ups = d)
for i_X = 2:N_X
    XiX  = Xs(i_X); upsd = upss(end);
    a  = (h_tau/(2*h_X^2))*upsd*XiX^2;
    td = (h_tau*kappa_tilde*(theta_tilde - upsd))/h_ups;
    e  = 1 - (h_tau/h_X^2)*upsd*XiX^2 + td;
    H_new(i_X, end) = a*(H(i_X-1, end) + H(i_X+1, end)) - td*H(i_X, end-1) + e*H(i_X, end);
end
end

function F = interpolacion_edp(S_pts, ups_pts, Xs, upss, H, r, T, b_dom, c_dom, d_dom)
S_pts   = S_pts(:);
ups_pts = ups_pts(:);
X_pts   = exp(r*T)*S_pts;
Xc      = min(max(X_pts,   0),           b_dom - 1e-9);
upsc    = min(max(ups_pts, c_dom + 1e-9), d_dom - 1e-9);
F       = zeros(length(S_pts), 1);
for k = 1:length(S_pts)
    F(k) = interp2(upss, Xs, H, upsc(k), Xc(k), 'linear')*exp(-r*T);
end
end

function [Fc, Fp] = ede_lote(S_pts, ups_pts, M_mc, n_mc, dt_mc, descuento, ...
                              K, kappa_tilde, theta_tilde, sigma, rho, r)
N_pts = length(S_pts);
blk   = 50;
n_b   = ceil(N_pts/blk);
Fc    = zeros(N_pts, 1);
Fp    = zeros(N_pts, 1);

for b = 1:n_b
    i0 = (b-1)*blk + 1;
    i1 = min(b*blk, N_pts);
    np = i1 - i0 + 1;
    Sb = S_pts(i0:i1).';
    Ub = ups_pts(i0:i1).';
    payoff_c = zeros(M_mc, np);
    payoff_p = zeros(M_mc, np);
    
    for m = 1:M_mc
        Sm = Sb; Um = Ub;
        for s = 1:n_mc
            dW1   = sqrt(dt_mc)*randn(1, np);
            dW2   = sqrt(dt_mc)*randn(1, np);
            ups_p = max(Um, 0);
            Sm    = Sm + r*Sm.*dt_mc + sqrt(ups_p).*Sm.*dW1;
            Um    = ups_p + kappa_tilde*(theta_tilde - ups_p)*dt_mc ...
                  + sigma*rho*sqrt(ups_p).*dW1 ...
                  + sigma*sqrt((1 - rho^2)*ups_p).*dW2;
        end
        payoff_c(m, :) = max(Sm - K, 0);
        payoff_p(m, :) = max(K - Sm, 0);
    end
    Fc(i0:i1) = descuento*mean(payoff_c, 1).';
    Fp(i0:i1) = descuento*mean(payoff_p, 1).';
end
end
