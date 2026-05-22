%%% Valoracion de Opciones Europeas en el Modelo de Heston
%%% Metodo: Semi--Analítico con Integración Numérica

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
divd    = 0;

% Parametros de la Opcion
K       = 80;
T       = 1;
S0      = 100;
ups0    = 0.04;

% Parametros del Modelo (bajo Q)
kappa_tilde = kappa + lambda*sigma;
theta_tilde = (kappa*theta - rho*sigma*(mu - r))/kappa_tilde;

% Fecha de Valoracion y Vencimiento
fecha_valo  = datenum(2026, 1 ,1);
fecha_venc  = datemnth(fecha_valo, round(T*12));

%% Precios de Call y Put

call = optByHestonNI(r, S0, fecha_valo, fecha_venc, 'call', K, ...
                    ups0, theta_tilde, kappa_tilde, sigma, rho, ...
                    'DividendYield', divd);
put = optByHestonNI(r, S0, fecha_valo, fecha_venc, 'put', K, ...
                    ups0, theta_tilde, kappa_tilde, sigma, rho, ...
                    'DividendYield', divd);

paridad_put_call = call - put - S0 + K*exp(-r*T);

tiempo = toc;

%% Superficie de Precios Call

b = 240; c = 0.01; d = 1;

n_S = 60; 
n_ups = 50;

Ss    = linspace(1, b*exp(-r*T), n_S);
upss  = linspace(c, d, n_ups);

for j = 1:n_ups
    C(:,j) = optByHestonNI(r, Ss, fecha_valo, fecha_venc, 'call', K, ...
                    upss(j), theta_tilde, kappa_tilde, sigma, rho, ...
                    'DividendYield', divd, 'ExpandOutput', true);

    P(:,j) = optByHestonNI(r, Ss, fecha_valo, fecha_venc, 'put', K, ...
                    upss(j), theta_tilde, kappa_tilde, sigma, rho, ...
                    'DividendYield', divd, 'ExpandOutput', true);
end

[ups_malla, S_malla] = meshgrid(upss, Ss);

figure(1); sub1 = subplot(1,2,1); hold on;
surf(ups_malla, S_malla, C); 
colormap(sub1, autumn); shading faceted; 
plot3(ups0, S0, call, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); 
xlabel('$\upsilon_t$', 'Interpreter', 'latex');
ylabel('$S_t$', 'Interpreter', 'latex');
zlabel('$C_t$', 'Interpreter', 'latex');
%title('Superficie de Precios Call');
%colorbar
view(-45, 30); grid on; hold off;

figure(1); sub2 = subplot(1,2,2); hold on;
surf(ups_malla, S_malla, P); 
colormap(sub2, summer); shading faceted; 
plot3(ups0, S0, put, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); 
xlabel('$\upsilon_t$', 'Interpreter', 'latex');
ylabel('$S_t$', 'Interpreter', 'latex');
zlabel('$P_t$', 'Interpreter', 'latex');
%title('Superficie de Precios Put');
%colorbar
view(-135, 30); grid on; hold off;

%% Mostrar Resultados

fprintf('==============================================================\n');
fprintf('    Valoracion de opciones europeas en el Modelo de Heston    \n');
fprintf('   Metodo en EDEs: Semi--Analítico con Integración Numérica   \n');
fprintf('==============================================================\n');
fprintf('Precio Call: %.6f.\n', call);
fprintf('Precio Put: %.6f.\n', put);
fprintf('Verificacion, Paridad Put-Call: %d.\n', paridad_put_call);
fprintf('Tiempo de Cómputo: %d.\n', tiempo);
fprintf('==============================================================\n');
