clc; clear all; close all;

%% Reducir muestras de 100 dimensiones a 2D
rng('shuffle');


% Cargo los patrones de muestra
datos = cell2mat(struct2cell(load("datos_para_clustering.mat")));

% Parametros
dim_espacio_entrada = size(datos, 2); 
dim_espacio_neuronas = 2;

p = size(datos, 1); % #patrones
xp = datos; % patrones

% Obtengo la muestra que tiene norma maxima.
norms = zeros(1, p);
for mu = 1:p
    norms(mu) = norm(datos(mu, :));
end
radio = max(norms);

n_x = 14; n_y = 14; % largo del eje x e y del espacio de coordenadas de neuronas.
n = n_x ^ dim_espacio_neuronas; % #neuronas

eta = 0.1; 
sigma = 8; % varianza inicial
alfa = 0.95; % cte. de actualizacion de la varianza
iteraciones = 50; % iteraciones maximas
iter = 1;


% inicializo los pesos
w = cell(sqrt(n), sqrt(n)); % vectores de pesos sinapticos de cada neurona. 

for i = 1:n_x
    for j = 1:n_y
        % Vectores columna de pesos, asociados a cada neurona.
        w{i, j} = rand(1, dim_espacio_entrada); 
        w{i, j} = radio .* w{i, j} ./ norm(w{i, j});
    end
end

%% 

% Entrenamiento de la red


neurona_ganadora = zeros(p, 2); % coordenadas de la neurona ganadora
dist_min = 1000 * ones(p, 1);
vecindad = zeros(p, n); % una matriz donde sus filas son listas con los valores de la vecindad de cada neurona; 
                        % cada fila es para un patron distinto

while(iter <= iteraciones)
    display("Iteración " + iter);
    for mu = shuffle(1:p)

        % 1. guardo el patron de entrada
        x = xp(mu, :)';
        display("patron "+mu);
        % 2,3. calculo la distancia; encuentro las neuronas ganadoras;
        dist = 0;
        dist_min = 1000 * ones(p, 1);
        for i = 1:n_x
            for j = 1:n_y

                dist = norm(x - w{i, j});

                if(dist <= dist_min(mu))
                    dist_min(mu) = dist;
                    neurona_ganadora(mu, :) = [i, j];
                end

            end
        end

        % 4,5. calculo las vecindades y actualizo los pesos
        idx_neu = 1;
        for i=1:n_x
            for j=1:n_y
                vecindad(mu, idx_neu) = func_vecindad([i, j], neurona_ganadora(mu, :), sigma); 
                w{i,j} = w{i,j} + eta * vecindad(mu, idx_neu) * (x - w{i,j});
                idx_neu = idx_neu + 1;
            end
        end


    end
    
    if(iter == 90 || iter == 50 || iter == 120 || iter == 140)

    end
    
    iter = iter + 1;
    
    sigma = sigma * alfa;
end

%% Armado de la función U



%%
plot_neuronas = figure();
vecindad_mu = reshape(vecindad(2, :), n_x, n_y);
imagesc(vecindad_mu);
title("Función de vecindad para el patrón " + mu , "interpreter", "tex");
hold on;
for i = 1:(n_x)
        plot([1, n_x], [i, i] , 'linewidth', 1, 'color', 'black' )   
        plot([i, i], [1, n_y] , 'linewidth', 1, 'color', 'black' )   
end
[X, Y] = meshgrid(1:n_x, 1:n_y);
scatter(X(:), Y(:), 100 , 'black', 'fill');

set(gca,'fontsize', 12);
set(plot_neuronas,'PaperSize', [20 10]); %set the paper size to what you want  
print(plot_neuronas,'Resultados/Ejercicio1/neuronas','-dpdf') % then print it

