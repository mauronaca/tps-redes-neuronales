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

n_x = 23; n_y = 23; % largo del eje x e y del espacio de coordenadas de neuronas.
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

%% Entrenamiento de la red 

% Entrenamiento de la red tarda como 2 horas! :O


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
                w{i,j} = w{i,j} + eta * vecindad(mu, idx_neu) * (x' - w{i,j});
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

w = struct2cell(load("w1.mat")); w= w{1}; % pesos de dimension 100!!!
w_mat = zeros(n_x, n_y); % matriz auxiliar que mapea las coordenadas con la matriz w de cell arrays.

for i = 1:n_x
 % workaround: en la linea de actualizacion de los pesos estaba haciendo (x
% - w{i,j}) y como x es un vector columna y w un vector fila, esa resta
% termina dando una matriz donde en cada columna esta el resultado replicado. El fix sería hacer (x'-w{i,j}), pero para los w
% que me guarde me quedo con los vector columna.

    for j = 1:n_y
        w{i,j} = w{i,j}(:,1);
    end
end
s = size(w);
B = zeros(s);
nn = numel(w);
matriz_vecinos = cell(s);

for ii = 1:nn
  B(ii) = 1;
  matriz_vecinos{ii} = w(bwdist(B,'ch') == 1);
  B(ii) = 0;
end

% Construyo la funcion U:
U = zeros(n_x, n_y); % matriz de la función U

for i = 1:n_x
    for j = 1:n_y
        
        % para la neurona [i,j]
        
        suma = 0; % inicializo la suma en 0 para luego guardarla en la matriz U.
        
        % recorro todos los vecinos, calculo la norma, y lo sumo en $suma.
        for vecinos = matriz_vecinos{i,j}
            for vecino_idx = 1:size(vecinos, 1)
                suma = suma + norm(w{i, j} - cell2mat(vecinos(vecino_idx)));
            end
        end
        
        U(i, j) = suma;
    end
end


% Grafico la matriz U

plot_clusters = figure();
pcolor(U);
colorbar;
set(gca,'fontsize', 12);
set(plot_clusters,'PaperSize', [20 10]); 
print(plot_clusters,'Resultados/Ejercicio3/clusters_1','-dpdf') 


