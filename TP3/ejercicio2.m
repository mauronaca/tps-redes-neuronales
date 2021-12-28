clc; clear all; close all;
 
%% Problema del viajante con red neuronal Kohonen.

rng('shuffle');

% Defino los patrones de entrada, muestras uniformes restringidas a un
% circulo de radio R. x(1)^2 + x(2)^2 -R = 0
p = 200; % #patrones
radio = 1; % Radio del circulo que contiene a los patrones
N = 2; % dim. de entrada
% genero patrones
%angulo = pi * linspace(0, 2, p); angulo = reshape(angulo, [p, 1]);
angulo = pi * unifrnd(0, 2, p, 1);
modulo = sqrt(unifrnd(0, radio, p, 1));
xp = [modulo .* cos(angulo), modulo .* sin(angulo)];
patrones_plot = figure();
viscircles([0,0], radio, 'color', 'black');
hold on;
scatter(xp(:, 1), xp(:, 2), 80, 'fill');

% Defino las neuronas y sus pesos. Los pesos se inicializan sobre el perimetro del circulo y de manera creciente
% por su angulo.
dim_espacio_neuronas = 1;
n_x = round(2.2 * p); n_y = 1; % largo del eje x e y del espacio de coordenadas de neuronas.
n = n_x ^ dim_espacio_neuronas; % #neuronas
w = cell(n_x); % vectores de pesos sinapticos de cada neurona. 
% inicializo los pesos
angulo = pi * linspace(0, 2, n);
for i = 1:n_x
    %for j = 1:n_y
    % Vectores columna de pesos, asociados a cada neurona.
    %angulo = pi * unifrnd(0, 2);
    modulo = sqrt(unifrnd(0, radio));
    w{i} = radio .* [ cos(angulo(i)) ; sin(angulo(i)) ]; 
    % grafico pesos
    scatter(w{i}(1), w{i}(2), 'fill', 'r');
    %end
end


set(gca,'fontsize', 12);
title(p + " ciudades y " + n + " neuronas");
xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
set(patrones_plot,'PaperSize', [20 10]); %set the paper size to what you want  
print(patrones_plot,'Resultados/Ejercicio2/patrones','-dpdf') % then print it

% Entreno la red

eta = 0.1;
sigma_init = 8;
sigma = sigma_init;
alfa = 0.95; % constante que decrementa la varianza

%iteraciones = size(sigma, 2);
iteraciones = 100;

neurona_ganadora = zeros(p, 2); % coordenadas de la neurona ganadora
vecindad = zeros(p, n); % una matriz donde sus filas son listas con los valores de la vecindad de cada neurona; 
                        % cada fila es para un patron distinto

iter = 1;


while(iter <= iteraciones)
    display("Iteración: " + iter);
    for mu = shuffle(1:p)

        % 1. guardo el patron de entrada
        x = xp(mu, :)';

        % 2,3. calculo la distancia; encuentro las neuronas ganadoras;
        dist = 0;
        dist_min = 1000 * ones(p, 1);
        for i = 1:n_x
            

            dist = norm(x - w{i});

            if(dist <= dist_min(mu))
                dist_min(mu) = dist;
                neurona_ganadora(mu, :) = i;
            end


        end

        % 4,5. calculo las vecindades y actualizo los pesos
        idx_neu = 1;
        for i=1:n_x

            vecindad(mu, idx_neu) = func_vecindad(i, neurona_ganadora(mu, :), sigma); 
            w{i} = w{i} + eta * vecindad(mu, idx_neu) * (x - w{i});
            idx_neu = idx_neu + 1;
            
        end


    end
    
    % Grafico la evolución de los pesos
    
    if(iter == 1 || iter == 10 || iter == 50)
        plot_ciudades = figure();
        hold on;
        viscircles([0,0], radio , 'color', 'black', 'linewidth', 0.5);
        scatter(xp(:, 1), xp(:, 2), 80 ,'black' );

        for i = 1:n_x

            % Vectores columna de pesos, asociados a cada neurona.
            % grafico pesos

            if i ~= n_x
                plot( [w{i + 1}(1), w{i}(1)], [w{i + 1}(2), w{i}(2)] ,'r' , 'linewidth', 1.5, 'color', 'black' )
            end
            scatter(w{i}(1), w{i}(2), 20 , 'fill', 'red');

        end
        % Uno el peso de la primera neurona con la ultima
        plot([w{1}(1), w{n}(1)], [w{1}(2), w{n}(2)], 'r', 'linewidth', 1, 'color', 'black');
        plot_ciudades.Renderer='Painters';

        title("\eta = " + eta + ", #neuronas = " + n + ", iteración = " + iter, 'interpreter', 'tex')
        set(gca,'fontsize', 12);
        xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
        ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
        set(plot_ciudades,'PaperSize', [20 10]); %set the paper size to what you want  
        path = "Resultados/Ejercicio2/recorrido" + iter;
        print(plot_ciudades,path,'-dpdf') % then print it


    end
    
    iter = iter + 1;
    
    % Actualizo el ancho de la vecindad
    sigma = sigma * alfa;
    
end

%%
plot_ciudades = figure();

hold on;

viscircles([0,0], radio , 'color', 'black', 'linewidth', 0.5);
scatter(xp(:, 1), xp(:, 2), 80 , 'fill' ,'black' );

for i = 1:n_x

    % Vectores columna de pesos, asociados a cada neurona.
    % grafico pesos

    if i ~= n_x
        plot( [w{i + 1}(1), w{i}(1)], [w{i + 1}(2), w{i}(2)] ,'r' , 'linewidth', 1.5, 'color', 'black' )
    end
    scatter(w{i}(1), w{i}(2), 20 , 'fill', 'red');

end
% Uno el peso de la primera neurona con la ultima
plot([w{1}(1), w{n}(1)], [w{1}(2), w{n}(2)], 'r', 'linewidth', 1, 'color', 'black');
plot_ciudades.Renderer='Painters';

title("\eta = " + eta + ", #neuronas = " + n + ", iteración = " + iter, 'interpreter', 'tex')
set(gca,'fontsize', 12);
xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
set(plot_ciudades,'PaperSize', [20 10]); %set the paper size to what you want  
print(plot_ciudades,'Resultados/Ejercicio2/recorrido','-dpdf') % then print it



plot_neuronas = figure();
hold on;
plot([1, n], [0, 0], 'black')
[X, Y] = meshgrid(1:n_x, 1:n_y);
scatter(1:1:n ,zeros(n,1), 100 , 'black', 'fill');
xlim([1, n]);
xlabel("Neuronas");
set(gca,'ytick',[]); % 1D
set(gca,'fontsize', 12);
set(plot_neuronas,'PaperSize', [20 10]); %set the paper size to what you want  
print(plot_neuronas,'Resultados/Ejercicio1/neuronas_rect','-dpdf') % then print it

