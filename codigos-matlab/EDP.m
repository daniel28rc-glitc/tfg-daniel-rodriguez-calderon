%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Metodo en EDP: Euler Explicito con Diferencias Finitas 

close all; clear all; clc; tic;

%% Parametros del Modelo

% Parametros del Modelo (bajo P)
r       = 0.05;
kappa   = 2;
theta   = 0.04;
sigma   = 0.3;
rho     = -0.5;
mu      = r;
lambda  = 0;

% Parametros de la Opcion
K       = 80;
T       = 1;
S0      = 100;
ups0    = 0.04;

% Parametros del Modelo (bajo Q)
kappa_tilde = kappa + lambda*sigma;
theta_tilde = (kappa*theta - rho*sigma*(mu - r))/kappa_tilde;

%% Dominio Truncado y Mallas Uniformes

b       = 240;
c       = 0.01;
d       = 1;
tau_max = T;

N_X     = 200;
N_ups   = 100;
N_tau   = 100000;

h_X     = b/N_X;
h_ups   = (d - c)/N_ups;
h_tau    = tau_max/N_tau;

X_vec = linspace(0, b, (N_X + 1)).';
ups_vec = linspace(c, d, (N_ups + 1));

%% Condiciones Iniciales y de Frontera

X0 = exp(r*T)*S0;

payoff_call = max((X_vec - K), 0);
payoff_put  = max((K - X_vec), 0);

H_call  = repmat(payoff_call, 1, (N_ups + 1));
H_put   = repmat(payoff_put, 1, (N_ups + 1));

H_call(1, :)    = 0;
H_call(end, :)  = max((b - K), 0);
H_put(1, :)     = K;
H_put(end, :)   = max((K - b), 0);

%% Aplicacion Euler Explicito e Interpolación del Precio

for i_tau = 1:N_tau
    H_call = paso_euler_edp(H_call, X_vec, ups_vec, h_tau, h_X, h_ups, ...
                             sigma, kappa_tilde, theta_tilde, rho);
    H_put = paso_euler_edp(H_put, X_vec, ups_vec, h_tau, h_X, h_ups, ...
                             sigma, kappa_tilde, theta_tilde, rho);
end

V_call = exp(-r*T)*H_call;
V_put = exp(-r*T)*H_put;

precio_call = interp_precio_2D(H_call, X_vec, ups_vec, X0, ups0, r, T);
precio_put  = interp_precio_2D(H_put, X_vec, ups_vec, X0, ups0, r, T);

paridad_put_call = precio_call - precio_put - S0 + K*exp(-r*T);

tiempo = toc;

%% Creacion de Figuras
[ups_mat, X_mat] = meshgrid(ups_vec, X_vec);
S_mat = exp(-r*T)*X_mat;

figure(1); sub1 = subplot(1,2,1); hold on;
surf(ups_mat, S_mat, V_call); colormap(sub1, autumn);  shading faceted;
plot3(ups0, S0, precio_call, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); 
xlabel('$\upsilon$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$S$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$C^{\textrm{EDP}}(0,S,\upsilon)$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18); 
%colorbar;
%title('Superficie de Precios Call');
view(-45, 30); grid on; hold off;

figure(1); sub2 = subplot(1,2,2); hold on;
surf(ups_mat, S_mat, V_put); colormap(sub2, summer); shading faceted; 
plot3(ups0, S0, precio_put, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); 
xlabel('$\upsilon$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$S$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$P^{\hspace{0.3mm} \textrm{EDP}}(0,S,\upsilon)$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 18); 
%colorbar;
%title('Superficie de Precios Put');
view(-135, 30); grid on; hold off;

%% Mostrar Resultados

fprintf('==============================================================\n');
fprintf('    Valoracion de opciones europeas en el Modelo de Heston    \n');
fprintf('    Metodo en EDP: Euler Explicito con Diferencias Finitas    \n');
fprintf('==============================================================\n');
fprintf('Precio Call: %.6f.\n', precio_call);
fprintf('Precio Put: %.6f.\n', precio_put);
fprintf('Verificacion, Paridad Put-Call: %d.\n', round(abs(paridad_put_call)));
fprintf('Tiempo de Cómputo: %d.\n', tiempo);
fprintf('==============================================================\n');

%% Funciones Locales

function H_new = paso_euler_edp(H, X_vec, ups_vec, h_tau, h_X, h_ups, ...
                        sigma, kappa_tilde, theta_tilde, rho)

    % Ajuste de Indices
    N_X = length(X_vec) - 1;
    N_ups = length(ups_vec) -1;
    H_new = H;

    % Nodos Interiores
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
            e_iX_iups = 1 - (h_tau/(h_X^2))*upsiups*XiX^2 ...
                        - (h_tau*sigma^2/(h_ups^2))*upsiups;

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
        XiX = X_vec(i_X);
        upsc = ups_vec(1);

        a_iX_cups = (h_tau/(2*h_X^2))*upsc*XiX^2;
        aux3      = (h_tau*kappa_tilde*(theta_tilde - upsc))/(h_ups);
        e_iX_cups = 1 - (h_tau/(h_X^2))*upsc*XiX^2 - aux3;

        H_new(i_X, 1) = a_iX_cups*(H(i_X-1, 1) + H(i_X+1, 1)) ...
                            + aux3*H(i_X,2) ...
                            + e_iX_cups*H(i_X, 1);
        
    end

    % Columna i_ups = N_ups + 1 (ups = d)
    for i_X = 2:N_X
        XiX = X_vec(i_X);
        upsd = ups_vec(end);

        a_iX_dups = (h_tau/(2*h_X^2))*upsd*XiX^2;
        aux3      = (h_tau*kappa_tilde*(theta_tilde - upsd))/(h_ups);
        e_iX_dups = 1 - (h_tau/(h_X^2))*upsd*XiX^2 + aux3;

        H_new(i_X, end) = a_iX_dups*(H(i_X-1, end) + H(i_X+1, end)) ...
                            - aux3*H(i_X,end-1) ...
                            + e_iX_dups*H(i_X, end);
        
    end
end

function precio = interp_precio_2D(H, X_vec, ups_vec, X0, ups0, r, T)
    iX0 = find(X_vec <= X0, 1, 'last');
    if isempty(iX0) || iX0 >= length(X_vec)
        iX0 = length(X_vec) - 1;
    end
    iX1 = iX0 + 1;
    tX  = (X0 - X_vec(iX0)) / (X_vec(iX1) - X_vec(iX0));

    iups0 = find(ups_vec <= ups0, 1, 'last');
    if isempty(iups0) || iups0 >= length(ups_vec)
        iups0 = length(ups_vec) - 1;
    end
    iups1 = iups0 + 1;
    tups = (ups0 - ups_vec(iups0)) / (ups_vec(iups1) - ups_vec(iups0));

    H_int = (1 - tX)*(1 - tups)*H(iX0, iups0) ...
            + tX*(1 - tups)*H(iX1, iups0)...
            + (1 - tX)*tups*H(iX0, iups1) + tX*tups*H(iX1, iups1);

    precio = exp(-r*T)*H_int;
end