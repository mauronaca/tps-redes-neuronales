clear all; clc; close all;
rng('shuffle');

%% Implementar RN que aprenda una funcion f(x,y,z)
% Correr esta seccion del script antes de correr alguna de las siuguientes.
rng('shuffle');

% Red neuronal de 2 capas, 16 neuronas en la capa oculta y 1 neurona en la
% capa de salida. Entrada de tamanio 3.
% 
%
%

% Set entrenamiento y test
p_train = 2000; % #patrones de entrenamiento
p_test = 1000;
p = p_train + p_test;
% Cada fila es un patron
X = [ones(p, 1), linspace(0, 2*pi, p)', linspace(0, 2*pi, p)', linspace(-1, 1, p)']; 
X(:, 2) = shuffle(X(:, 2)); X(:, 3) = shuffle(X(:, 3)); X(:, 4) = shuffle(X(:, 4));
X_train = X(1:p_train, :); X_test = X(1:p_test, :);

% Salida
Yd = sin(X(:, 2)) + cos(X(:, 3)) + X(:, 4);
Yd_train = sin(X_train(:, 2)) + cos(X_train(:, 3)) + X_train(:, 4);
Yd_test = sin(X_test(:, 2)) + cos(X_test(:, 3)) + X_test(:, 4);

% Normalizo a -1 y 1:
normY = 6;
Yd_train = Yd_train / normY;
Yd_test = Yd_test / normY;
normX = 2*pi;
%X_train(:,2:4) = X_train(:,2:4) / normX;
%X_test(:,2:4) = X_test(:,2:4) / normX;

Y_train = zeros(p_train, 1); 
Y_test = zeros(p_test, 1);

M = 2; % #capas
n = 3 + 1; % #entradas
N = [16, 1]; % #neuronas p/capa


% Inicializo los pesos:
w = { randi([-10, 10], N(1), n)/1000, randi([-10, 10], N(2), N(1))/1000 };
dw = { zeros(N(1), n), zeros(N(2), N(1)) };

% bias
b = {0, 1};

% Estados de las neuronas:
h = {zeros(N(1), 1), zeros(N(2), 1)};
V = {X_train(1, :)', tanh(h{1}), tanh(h{2})};
% Error de cada capa
sigma = {zeros(N(1), 1), zeros(N(2), 1)};

tolerancia = 0.001;
epochMax = 1000;
error_train = 1;
errors_train = zeros(1, epochMax); % Guardo el error en cada epoca para luego graficarlo en funcion de las epocas.
errors_test = zeros(1, epochMax); 
error_test = 1;
epoch = 0;
iter = 0;
mini_batch_size = 1;
cant_batch = (p_train / mini_batch_size);

constAprendizaje = 0.01;

%% Entreno y grafico error vs epoca
while (error_test > tolerancia) && (epoch < epochMax)
    
    epoch = epoch + 1;
    %iter = 0;
    % Pongo en cero el delta w.
    %dw = { zeros(N(1), n), zeros(N(2), N(1)) };
    
    % Itero sobre todos los mini batches:     
    for batch = randperm(cant_batch)
        %disp("Batch: " + batch)
        dw = { zeros(N(1), n), zeros(N(2), N(1)) };
        for mu = batch:(batch + mini_batch_size - 1)
            % Itero sobre los patrones del mini batch. Para el punto (a) es
            % 1 solo patron.
            % 1- Calculo la salida real ante el patrón de entrada. 2-
            % Calculo error y hago backrpoagation. 3- Calculo los dw y los
            % acumulo. 
            
            x = X_train(mu, :)'; % Selecciono el patrón.
            V{1} = x;
            
            % Calculo la salida real
            for m = 1:M
                % Itero sobre todas las capas
                h{m} = w{m} * V{m} + b{m};
                V{m+1} = tanh(h{m});                
            end
            
            % Backpropagation: Calculo el error desde la capa de salida hacia
            % la de entrada:
            % Primero se calcula el error entre la salida y la entrada en la
            % ultima capa:
            sigma{M} = (1 - tanh(h{M})^2) * (Yd_train(mu) - V{M+1});
            % las demas capas:
            for m  = M-1:-1:1
                sigma{m} = (1 - tanh(h{m}).^2) .* (w{m+1} * sigma{m+1})';
            end
            
            
            % Calculo los delta w
            for m=1:M
                dw{m} = dw{m} + constAprendizaje * sigma{m} * V{m}';
            end
            
        end
        
        % Actualizo w.
        for m = 1:M
            w{m} = w{m} + dw{m};
        end
        
        iter = iter + 1;
    end
    
    % Calculo el error de testeo usando los patrones de test.
    % Computo la salida
    for mu = 1:p_test
        Y_test(mu) = computarSalida(X_test(mu, :)', w, b);
    end
    for mu = 1:p_train
        Y_train(mu) = computarSalida(X_train(mu, :)', w, b);
    end
    
    error_test = sum((Yd_test - Y_test) .^ 2)/p_test;
    errors_test(epoch) = error_test;
    error_train = sum((Yd_train - Y_train) .^ 2)/p_train;
    errors_train(epoch) = error_train;
    
    disp("Epoch: " + epoch + "; Error de testeo: " + errors_test(epoch) + "; Error de training: " + errors_train(epoch));
end

if (epoch == epochMax)
    disp("No convergió, error " + errors_test(epochMax) + ", iteraciones: " + iter);
else
    disp("Cantidad de epocas " + epoch + " y error " + errors_test(epoch));
end

plot_error_epoch = figure();
plot(1:epochMax, errors_test, 'linewidth', 2.5, 'color', 'black');
legend("Tamaño mini batch: " + mini_batch_size);
xlabel('Epoca'); ylabel('Error');
set(gca,'fontsize', 12);
xlim([1 epoch])
set(plot_error_epoch, 'PaperSize',[20 10]); 
%print(plot_error_epoch,'Resultados/Ejercicio4/puntoA.pdf','-dpdf');

pause;

%% Punto B) Vario el tamaño de los batches, y la const de aprendizaje. Grafico el tiempo y epocas en f de eso.
rng('shuffle');

iter_vs_batch_size = zeros(1, size(batchSizes(100, 200), 2));
tiempo_vs_batch_size = zeros(1, size(batchSizes(100, 200), 2));
eta_idx = 1;
w_init = { randi([-10, 10], N(1), n)/1000, randi([-10, 10], N(2), N(1))/1000 };

for mini_batch_size = batchSizes(100, p_train)

    disp("Tamaño mini batch: " + mini_batch_size);
    error_test = 1;
    epoch = 0;
    w = w_init;
    iter = 0;
    tic;
    
    while (error_test > tolerancia) && (epoch < epochMax)

        epoch = epoch + 1;

        % Itero sobre todos los mini batches:     
        for batch = randperm((p_train / mini_batch_size))
            dw = { zeros(N(1), n), zeros(N(2), N(1)) };
            
            for mu = batch:(batch + mini_batch_size - 1)
                x = X_train(mu, :)'; 
                V{1} = x;
                for m = 1:M
                    h{m} = w{m} * V{m} + b{m};
                    V{m+1} = tanh(h{m});
                end
                sigma{M} = (1 - tanh(h{M})^2) * (Yd_train(mu) - V{M+1});
                for m  = M-1:-1:1
                    sigma{m} = (1 - tanh(h{m}).^2) .* (w{m+1} * sigma{m+1})';
                end
                for m=1:M
                    dw{m} = dw{m} + constAprendizaje * sigma{m} * V{m}';
                end
            end

            for m = 1:M
                w{m} = w{m} + dw{m};
            end

            iter = iter + 1;
        end

        for mu = 1:p_test
            Y_test(mu) = computarSalida(X_test(mu, :)', w, b);
        end

        error_test = sum((Yd_test - Y_test) .^ 2)/p_test;
        disp("Epoch: " + epoch + "; Error de testeo: " + error_test );
    end
    
    iter_vs_batch_size(eta_idx) = iter;
    tiempo_vs_batch_size(eta_idx) = toc;
    eta_idx = eta_idx + 1;
end

% Plot de iteraciones
plot_iter_batch = figure();
plot(batchSizes(100, p_train), iter_vs_batch_size, 'linewidth', 2.5, 'color', 'black');
xlabel('Tamaño batch'); ylabel('Iteraciones');
legend("eta = " + constAprendizaje);
set(gca,'fontsize', 12);
set(plot_iter_batch, 'PaperSize', [20 10]); 
%print(plot_iter_batch,'Resultados/Ejercicio4/iterVSbatch.pdf','-dpdf');

% Plot del tiempo
plot_tiempo_batch = figure();
plot(batchSizes(100, p_train), tiempo_vs_batch_size, 'linewidth', 2.5, 'color', 'black');
xlabel('Tamaño batch'); ylabel('Tiempo');
legend("eta = " + constAprendizaje);
set(gca,'fontsize', 12);
set(plot_tiempo_batch, 'PaperSize', [20 10]); 
%print(plot_tiempo_batch,'Resultados/Ejercicio4/tiempoVSbatch.pdf','-dpdf');

pause;
%% Punto b) Graficar iter vs constAprendizaje
rng('shuffle');

etas = [0.01:0.01:0.09, 0.1:0.1:0.9];
iter_vs_eta = zeros(1, size(etas, 2));
time_vs_eta = zeros(1, size(etas, 2));
eta_idx = 1;
w_init = { randi([-10, 10], N(1), n)/1000, randi([-10, 10], N(2), N(1))/1000 };
mini_batch_size = 1;


for eta = etas

    disp("Const. Aprendizaje: " + eta);
    error_test = 1;
    epoch = 0;
    w = w_init;
    iter = 0;
    tic;
    
    while (error_test > tolerancia) && (epoch < epochMax)

        epoch = epoch + 1;

        for batch = randperm((p_train / mini_batch_size))
            dw = { zeros(N(1), n), zeros(N(2), N(1)) };
            
            for mu = batch:(batch + mini_batch_size - 1)
                x = X_train(mu, :)'; 
                V{1} = x;
                for m = 1:M
                    h{m} = w{m} * V{m} + b{m};
                    V{m+1} = tanh(h{m});
                end
                sigma{M} = (1 - tanh(h{M})^2) * (Yd_train(mu) - V{M+1});
                for m  = M-1:-1:1
                    sigma{m} = (1 - tanh(h{m}).^2) .* (w{m+1} * sigma{m+1})';
                end
                for m=1:M
                    dw{m} = dw{m} + eta * sigma{m} * V{m}';
                end
            end

            for m = 1:M
                w{m} = w{m} + dw{m};
            end

            iter = iter + 1;
        end

        for mu = 1:p_test
            Y_test(mu) = computarSalida(X_test(mu, :)', w, b);
        end

        error_test = sum((Yd_test - Y_test) .^ 2)/p_test;
        disp("Epoch: " + epoch + "; Error de testeo: " + error_test );
    end
    
    iter_vs_eta(eta_idx) = iter;
    time_vs_eta(eta_idx) = toc;
    eta_idx = eta_idx + 1;
end


% Plot de iter vs eta
plot_eta_iter = figure();
plot(etas, iter_vs_eta, 'linewidth', 2.5, 'color', 'black');
xlabel('Constante de Aprendizaje'); ylabel('Iteraciones');
legend("Tamaño batch = " + mini_batch_size);
set(gca,'fontsize', 12);
set(plot_eta_iter, 'PaperSize', [20 10]); 
print(plot_eta_iter,'Resultados/Ejercicio4/iterVSeta.pdf','-dpdf');


% Plot de tiempo vs eta
plot_time_iter = figure();
plot(etas, time_vs_eta, 'linewidth', 2.5, 'color', 'black');
xlabel('Constante de Aprendizaje'); ylabel('Tiempo');
legend("Tamaño batch = " + mini_batch_size);
set(gca,'fontsize', 12);
set(plot_time_iter, 'PaperSize', [20 10]); 
print(plot_time_iter,'Resultados/Ejercicio4/timeVSeta.pdf','-dpdf');


%% Funciones utiles

function y = salida(x, w, b, normx, normy)
    x(2:4) = x(2:4) / normx;
    y = computarSalida(x, w, b) * normy;
end

function y = batchSizes(max, p)
    iter = 1;
    for i = 1:max
        if rem(p, i) == 0
            y(iter) = i;
            iter = iter + 1;
        end
    end
end

