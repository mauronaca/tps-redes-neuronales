clear all; clc; close all;
rng('shuffle');

%% Red multicapa que resuelve XOR de 2 entradas mediante Simulated Annealing

M = 2; % #capas
N = [3, 1]; % #neuronas
n = 2 + 1; % #entradas + bias
p_train = 4; % #patrones

% Patrones
X_train = [1, -1, -1; 1, -1, 1; 1, 1, -1; 1, 1, 1];
% Salida deseada y salida de la red
Yd_train = [-1; 1; 1; -1];
Y = zeros(p_train, 1);
Y_estrella = zeros(p_train, 1);


% Inicializo estado de las neuronas y salidas
h = {zeros(N(1), 1), zeros(N(2), 1)};
V = {X_train(1, :)', tanh(h{1}), tanh(h{2})};

% Inicializo pesos y sus deltas
w = { rand(N(1), n) - 1/2, ( rand(N(2), N(1)) - 1/2 )}; % Pesos para cada capa
w_estrella = { zeros(N(1), n), zeros(N(2), N(1)) };
dw = { zeros(N(1), n), zeros(N(2), N(1)) };

% bias
b = {(zeros(N(1), 1)), (rand(N(2), 1)-0.5)};

% Constantes
sigma = rand * 0.8;
beta = 1;
alfa = 0.8; 
T = 100; 

max_iter = 1000; % Corte por maximo de iteaciones
tolerancia = 0.01;

error = 1;
E = ones(1, max_iter);
E_estrella = ones(1, max_iter);
dE = ones(1, max_iter);

%% Entrenamiento

k = 0; % iteracion

while (error > tolerancia) && (k < max_iter)
    
    k = k + 1;
    disp("Iteracion: " + k);

    % Simulated Annealing:
    % 1. calculo la salida con w
    % 2. calculo el error sobre todos los patrones utilizando la salida
    % calculada anteriormente
    % 3. genero w* con muestras
    % 4. obtengo la salida con w* y computo E*
    % 5. actualizo w segun E-E_estrella
    % 6. actualizo T
    
    % 1. calculo la salida con w para todos los patrones
    for mu = 1:p_train
        X = X_train(mu, :)';
        Y(mu) = computarSalida(X, w, b);
    end
    
    % 2. calculo el E con w.
    E(k) = sum((Yd_train - Y) .^ 2) / p_train;
    error = E(k); % Guardo el error actual
    disp("Error: " + error);
    
    % 3. genero w*
    for m = 1:M
        for i = randperm(size(w{m}, 1))
            for j = randperm(size(w{m}, 2))
                dw{m}(i,j) = normrnd(0, (sigma));
                w_estrella{m}(i,j) = w{m}(i,j) + dw{m}(i,j);
            end
        end

    end
    
    % 4. obtengo el E*
    for mu = 1:p_train
        X = X_train(mu, :)';
        Y_estrella(mu) = computarSalida(X, w_estrella, b);
    end
    E_estrella(k) = sum((Yd_train - Y_estrella) .^ 2) / p_train;
    
    % 5. actualizo w
    % Compouto la diferencia entre los errores
    dE(k) = E_estrella(k) - E(k);
    p = exp(-dE(k) / (beta * T)); % Probabilidad
    r = rand;
    
    if(dE(k) <= 0)
       % Acepto el cambioo
       w = w_estrella; 
       E(k) = E_estrella(k);
       error = E_estrella(k);
    else
        if(r < p)
           w = w_estrella; 
           E(k) = E_estrella(k);
           error = E_estrella(k);
        end
    end
    
    % 6. actualizo T
    T = alfa * T ;
    %sigma = sigma * alfa;
    
end

%% Plot

grafico = figure();
plot(1:k, E(1:k), 'linewidth', 3, 'color', 'black');
xlabel('Iteraciones'); ylabel('Error');
titulo = "Simulated Annealing con $T_{inicial}$ = " + T + ", $\alpha$ = " + alfa + ", $\sigma$ = " + sigma;
title(titulo,'Interpreter','latex')
legend("T=100")
set(gca,'fontsize', 12);
grid('on');
set(grafico, 'PaperSize',[25 10]); %set the paper size to what you want  
print(grafico, 'Resultados/Ejercicio6/grafico6-1','-dpdf', '-bestfit') % then print it
