%%% SIMULACION DE TRAYECTORIAS DE MOVIMIENTOS BROWNIANOS
clc; clear all; close all;

% Parametros de la simulacion
T = 1;
N = 1000;
dt = T/N;
t = linspace(0,T,N+1);

% Numero de caminos en cada subplot
trays = [1, 5, 20, 100];
max_trays = max(trays);

% Color de representación
color_escuela = [0.0, 0.25, 0.52];

% Generar trayectorias a dibujar
rng(28); % Semilla fija para reproducibilidad
dW = sqrt(dt)*randn(max_trays, N);
B = zeros(max_trays, N+1);
B(:, 2:end) = cumsum(dW, 2);

% Realizar los plots
figura = figure('Position', [100, 100, 800, 800], 'Color', 'w');

% Bucle de representación de los subplots
for k = 1:4
    subplot(2, 2, k);
    hold on;

    % Determinar la opacidad y el grosor dependiendo de la cantidad 
    % de trayectorias
    if trays(k) <= 5
        opac = 1.0;
        gros = 1.0;
    elseif trays(k) == 20
        opac = 0.6;
        gros = 0.8;
    else
        opac = 0.3;
        gros = 0.5;
    end

    % Graficar las trayectorias
    for i = 1:trays(k)
        plot(t, B(i,:), 'Color', [color_escuela, opac], 'LineWidth', gros);
    end

    hold off

    % Estilo y etiquetas 
    box on; grid on;
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$B_t$', 'Interpreter', 'latex', 'FontSize', 11);
    xlim([0, 1]);
    if trays(k) == 1
        texto = sprintf('%d trayectoria', trays(k));
    else
        texto = sprintf('%d trayectorias', trays(k));
    end
    title(texto, ...
             'Units', 'normalized', 'Position', [0.5, 0.1, 0], ...
             'BackgroundColor', 'w', 'EdgeColor', 'k', ...
             'FontWeight', 'normal', 'FontSize', 10);
end

% Añadir el titulo general 
titulo = sgtitle('', ...
        'FontWeight', 'bold', 'FontSize', 15);
title.Color = color_escuela; 
