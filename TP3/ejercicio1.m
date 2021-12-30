clc; clear all; close all

%% Red Kohnonen de 2 entradas, distribución dentro de un circulo
rng('shuffle');

% Defino los patrones de entrada, muestras uniformes restringidas a un
% circulo de radio R. x(1)^2 + x(2)^2 -R = 0
p = 50; % #patrones
radio = 1; % Radio del circulo que contiene a los patrones
N = 2; % dim. de entrada
% genero patrones
modulo = sqrt(unifrnd(0, radio, p, 1));
angulo = pi * unifrnd(0, 2, p, 1);
L = 1; B = 2;
a1 = 0;
b1 = L;
a = (b1-a1).*rand(p,1) + a1;
a2 = 0;
b2 = B;
b = (b2-a2).*rand(p,1) + a2;
%xp = [a, b]; 
xp = [modulo .* cos(angulo), modulo .* sin(angulo)];
patrones_plot = figure();
viscircles([0,0], radio, 'color', 'black');
hold on;
Rect = [ 0 0 ; L 0 ; L B ; 0 B ; 0 0] ;
%plot(Rect(:,1),Rect(:,2), 'color', 'black')
%scatter(xp(:, 1), xp(:, 2), 'fill');
scatter(xp(:, 1), xp(:, 2), 'fill');

% Defino las neuronas y sus pesos
dim_espacio_neuronas = 2;
n_x = 14; n_y = 14; % largo del eje x e y del espacio de coordenadas de neuronas.
n = n_x ^ dim_espacio_neuronas; % #neuronas
w = cell(sqrt(n), sqrt(n)); % vectores de pesos sinapticos de cada neurona. 
% inicializo los pesos
for i = 1:n_x
    for j = 1:n_y
        % Vectores columna de pesos, asociados a cada neurona.
        angulo = pi * unifrnd(0, 2);
        modulo = sqrt(unifrnd(0, radio));
        w{i, j} = modulo .* [ cos(angulo) ; sin(angulo) ]; 
        % grafico pesos
        scatter(w{i,j}(1), w{i,j}(2), 'fill', 'r');

    end
end


set(gca,'fontsize', 12);
xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
set(patrones_plot,'PaperSize', [20 10]); %set the paper size to what you want  
print(patrones_plot,'Resultados/Ejercicio1/patrones','-dpdf') % then print it


% Pasos algoritmo
% 1. se elige un patron x al azar.
% 2. se calcula la distancia ese vector x y los vectores w de cada neurona.
% 3. obtengo la neurona ganadora, la cual tiene la menor distancia.
% 4. obtengo la funcion de vecindad 

%
eta = 0.1;
patron_pa_graficar = randi([1, p]);
%sigma = [150:-1:1];
sigma = 8;
alfa = 0.95;

iteraciones = 100;

neurona_ganadora = zeros(p, 2); % coordenadas de la neurona ganadora
dist_min = 1000 * ones(p, 1);
vecindad = zeros(p, n); % una matriz donde sus filas son listas con los valores de la vecindad de cada neurona; 
                        % cada fila es para un patron distinto

iter = 1;


while(iter <= iteraciones)
    display("Iteración " + iter);

    for mu = shuffle(1:p)

        % 1. guardo el patron de entrada
        x = xp(mu, :)';

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
    if(iter == 1 || iter == 10 || iter == 50)
        plot_pesos = figure();
        hold on;
        viscircles([0,0], radio , 'color', 'black', 'linewidth', 0.5);
        scatter(xp(:, 1), xp(:, 2), 50 ,'fill','blue' );

        for i = 1:n_x
            for j = 1:n_y
                % Vectores columna de pesos, asociados a cada neurona.
                % grafico pesos
                if j ~=	 n_y
                    plot( [w{i, j}(1), w{i, j+1}(1)], [w{i, j}(2), w{i, j+1}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
                end
                if i ~= n_x
                    plot( [w{i + 1, j}(1), w{i, j}(1)], [w{i + 1, j}(2), w{i, j}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
                end
                scatter(w{i, j}(1), w{i, j}(2), 50 , 'fill', 'red');
            end
        end
        % Uno el peso de la primera neurona con la ultima
        plot_pesos.Renderer='Painters';

        title("\eta = " + eta + ", #neuronas = " + n + ", iteración = " + iter, 'interpreter', 'tex')
        set(gca,'fontsize', 12);
        xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
        ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
        set(plot_pesos,'PaperSize', [20 10]); %set the paper size to what you want  
        path = "Resultados/Ejercicio1/pesos" + iter;
        print(plot_pesos,path,'-dpdf') % then print it


    end
    iter = iter + 1;
    sigma = sigma * alfa;

end

plot_pesos = figure();
plot_pesos.Renderer='Painters';
title("\eta = " + eta + ", #neuronas = " + n + ", #patrones = " + p + "iteracion " + iter, 'interpreter', 'tex')
viscircles([0,0], radio, 'color', 'black');
hold on;
scatter(xp(:, 1), xp(:, 2), 100 , 'fill' );

for i = 1:n_x
    for j = 1:n_y
        % Vectores columna de pesos, asociados a cada neurona.
        % grafico pesos
        if j ~=	 n_y
            plot( [w{i, j}(1), w{i, j+1}(1)], [w{i, j}(2), w{i, j+1}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
        end
        if i ~= n_x
            plot( [w{i + 1, j}(1), w{i, j}(1)], [w{i + 1, j}(2), w{i, j}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
        end
        scatter(w{i, j}(1), w{i, j}(2), 50 , 'fill', 'red');
    end
end

set(gca,'fontsize', 12);
xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
set(plot_pesos,'PaperSize', [20 10]); %set the paper size to what you want  
print(plot_pesos,'Resultados/Ejercicio1/pesos','-dpdf') % then print it



%% Red Kohnonen de 2 entradas, distribución dentro de un rectangulo
clc; clear all;
rng('shuffle');

% Defino los patrones de entrada, muestras uniformes restringidas a un
% circulo de radio R. x(1)^2 + x(2)^2 -R = 0
p = 50; % #patrones
N = 2; % dim. de entrada
% genero patrones
L = 1; B = 1;
a1 = 0;
b1 = L;
a = (b1-a1).*rand(p,1) + a1;
a2 = 0;
b2 = B;
b = (b2-a2).*rand(p,1) + a2;
xp = [a, b];
patrones_plot = figure();
hold on;
Rect = [ 0 0 ; L 0 ; L B ; 0 B ; 0 0] ;
plot(Rect(:,1),Rect(:,2), 'color', 'black','linewidth' , 2)
scatter(xp(:, 1), xp(:, 2), 'fill', 'b');

% Defino las neuronas y sus pesos
dim_espacio_neuronas = 2;
n_x = 14; n_y = 14; % largo del eje x e y del espacio de coordenadas de neuronas.
n = n_x ^ dim_espacio_neuronas; % #neuronas
w = cell(sqrt(n), sqrt(n)); % vectores de pesos sinapticos de cada neurona. 
% inicializo los pesos
for i = 1:n_x
    for j = 1:n_y
        % Vectores columna de pesos, asociados a cada neurona.
        a = (b1-a1).*rand + a1;
        b = (b2-a2).*rand + a2;
        w{i, j} = [ a ; b ]; 
        scatter(w{i,j}(1), w{i,j}(2), 'fill', 'r');
    end
end


set(gca,'fontsize', 12);
xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
set(patrones_plot,'PaperSize', [20 10]); %set the paper size to what you want  
print(patrones_plot,'Resultados/Ejercicio1/patrones_rect','-dpdf') % then print it


% Pasos algoritmo
% 1. se elige un patron x al azar.
% 2. se calcula la distancia ese vector x y los vectores w de cada neurona.
% 3. obtengo la neurona ganadora, la cual tiene la menor distancia.
% 4. obtengo la funcion de vecindad 

%
eta = 0.75;
patron_pa_graficar = randi([1, p]);
sigma = [200:-1:1];

iteraciones = size(sigma, 2);

neurona_ganadora = zeros(p, 2); % coordenadas de la neurona ganadora
dist_min = 1000 * ones(p, 1);
vecindad = zeros(p, n); % una matriz donde sus filas son listas con los valores de la vecindad de cada neurona; 
                        % cada fila es para un patron distinto

iter = 1;


while(iter <= iteraciones)
    display("Iteración " + iter);

    for mu = shuffle(1:p)

        % 1. guardo el patron de entrada
        x = xp(mu, :)';

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
                vecindad(mu, idx_neu) = func_vecindad([i, j], neurona_ganadora(mu, :), sigma(iter)); 
                w{i,j} = w{i,j} + eta * vecindad(mu, idx_neu) * (x - w{i,j});
                idx_neu = idx_neu + 1;
            end
        end


    end
    
    if(iter == 1 || iter == 10 || iter == 50 || iter == 150)
        plot_pesos = figure();
        hold on;
        plot(Rect(:, 1), Rect(:, 2), 'color', 'black', 'linewidth', 2)
        scatter(xp(:, 1), xp(:, 2), 50 ,'fill','blue' );

        for i = 1:n_x
            for j = 1:n_y
                % Vectores columna de pesos, asociados a cada neurona.
                % grafico pesos
                if j ~=	 n_y
                    plot( [w{i, j}(1), w{i, j+1}(1)], [w{i, j}(2), w{i, j+1}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
                end
                if i ~= n_x
                    plot( [w{i + 1, j}(1), w{i, j}(1)], [w{i + 1, j}(2), w{i, j}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
                end
                scatter(w{i, j}(1), w{i, j}(2), 50 , 'fill', 'red');
            end
        end
        % Uno el peso de la primera neurona con la ultima
        plot_pesos.Renderer='Painters';

        title("\eta = " + eta + ", #neuronas = " + n + ", iteración = " + iter, 'interpreter', 'tex')
        set(gca,'fontsize', 12);
        xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
        ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
        set(plot_pesos,'PaperSize', [20 10]); %set the paper size to what you want  
        path = "Resultados/Ejercicio1/pesos_rect" + iter;
        print(plot_pesos,path,'-dpdf') % then print it


    end
    
    iter = iter + 1;

end
%%
plot_pesos_rect = figure();

Rect = [ 0 0 ; L 0 ; L B ; 0 B ; 0 0] ;
plot(Rect(:, 1), Rect(:, 2), 'color', 'black', 'linewidth', 2)
hold on;
scatter(xp(:, 1), xp(:, 2), 50 , 'fill' , 'b');
title("\eta = " + eta + ", #neuronas = " + n + ", #patrones = " + p + " iteracion " + iter, 'interpreter', 'tex');

for i = 1:n_x
    for j = 1:n_y
        % Vectores columna de pesos, asociados a cada neurona.
        % grafico pesos
        if j ~=	 n_y
            plot( [w{i, j}(1), w{i, j+1}(1)], [w{i, j}(2), w{i, j+1}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
        end
        if i ~= n_x
            plot( [w{i + 1, j}(1), w{i, j}(1)], [w{i + 1, j}(2), w{i, j}(2)] ,'r' , 'linewidth', 1, 'color', 'black' )
        end
        scatter(w{i, j}(1), w{i, j}(2), 50 , 'fill', 'red');
    end
end

set(gca,'fontsize', 12);
xlabel('\color{blue}\xi(1) \color{red}w(1)', 'interpreter', 'tex');
ylabel('\color{blue}\xi(2) \color{red}w(2)', 'interpreter', 'tex');
set(plot_pesos_rect,'PaperSize', [20 10]); %set the paper size to what you want  
print(plot_pesos_rect, 'Resultados/Ejercicio1/pesos_rect','-dpdf') % then print it


%%


plot_neuronas = figure();


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

