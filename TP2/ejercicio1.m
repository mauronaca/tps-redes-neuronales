%------------------------------------------------------------------------------
% 1) Implemente un perceptrón simple que aprenda la función lógica AND de 2 y de 4
% entradas. Lo mismo para la función lógica OR. Para el caso de 2 dimensiones, grafique
% la recta discriminadora y todos los vectores de entrada de la red.
% ------------------------------------------------------------------------------
close all; clear all; clc;
rng('shuffle');

%% Perceptron simple para funcion AND de 2 entradas
% ------------------------------------------------
%    x1 ? x2  ? yd
%   ----------------
%   -1  ? -1  ? -1
%   -1  ?  1  ? -1
%    1  ? -1  ? -1
%    1  ?  1  ?  1

n = 3; % Cantidad de entradas
p = 4; % Canntidad de patrones


% Patrones. tamaño (4, 2) 4 patrones, 3 entradas
xp = [1, -1, -1; 1, -1, 1; 1, 1, -1; 1, 1, 1]; % Entrada 0 es 1
%xp = transpose(xp);
% Salida deseada (clasificacion 1 o -1)
yd = [-1; -1; -1; 1];

% Inicio w al azar. Como son 3 entradas w es de tamaño (3,1). w =  w0,w1, w2
%w = rand(n, 1);
w = zeros(n, 1);
error = 1;
y = zeros(p, 1);
cantIter = 0;
constAprendizaje = 0.9;

while (error > 0)
    
    % Elijo un patron al azar
    mu = randi(4);
    disp("Patron " + mu);
    
    % Me quedo con dicho patron y se lo asigno a x
    x = transpose(xp(mu, :));
    disp("Patron actual: " + x);
    disp("w actual: " + w);
    
    % Calculo la salida del perceptron ante la entrada x
    y(mu) = signo(transpose(w) * x);
    
    % Actaulizacion de w, aprendizaje.
    w = w + constAprendizaje * (yd(mu) - y(mu)) * x;
    
    error = computarError(y, yd);
    disp("Error " + error);
    
    cantIter = cantIter + 1;
    
end

w0 = w(1); w1 = w(2); w2 = w(3);
x1 = -2:0.1:2;
x2 = -w0 / w2 - x1 .* (w1 / w2);

display("Iteraciones: " + cantIter);
display(w);
display("y: " + y); display("yd: " + yd);

graficoAnd = figure();
plot(x1, x2, 'linewidth', 3, 'color', 'red')
%fplot(@(t) -w0/w2 - t*w1/w2);
hold on;

scatter(xp(:, 2), xp(:, 3), 80, 'filled', 'black')
xlim([-2, 2]); ylim([-2, 2]); 
gridxy(0,0, 'Linestyle',':'); 
gridxy(1,1, 'Linestyle',':'); 
gridxy(-1,-1,'Linestyle',':');
title("Perceptron Simple AND 2 entradas");
set(gca,'fontsize', 14);
set(graficoAnd,'PaperSize',[20 10]); %set the paper size to what you want  
%print(graficoAnd,'Resultados/Ejercicio1/and2entradas','-dpdf') % then print it

%% Saco un promedio de iteraciones en funcion de la const de Aprendizaje
for a = 1 % Para ocultar el codigo
    constantes = 0.1:0.05:1;
    cantIter = 200;
    tipoDeInicializacion = [rand(n, 1), zeros(n, 1), ones(n, 1)]; 
    cantInit = size(tipoDeInicializacion); cantInit = cantInit(2);
    legendTipoDeInit = [...
        "Inicializacion de w distribucion uniforme",... 
        "Inicializacion de w con 0s",...
        "Inicializacion de w con 1s"];

    grafico2 = figure();
    for init = 1:cantInit

        iteraciones = zeros(size(constantes));
        pasos = zeros(1, cantIter);
        iterConst = 1;

        for constAprendizaje = constantes

            for i = 1:cantIter

                % Inicio w al azar. Como son 3 entradas w es de tamaño (3,1). w =  w0,w1, w2
                w = tipoDeInicializacion(:, cantInit);
                error = 1;
                y = zeros(p, 1);
                cantPasos = 0;

                % Entreno el perceptron
                while (error > 0)

                    % Elijo un patron al azar
                    mu = randi(4);

                    % Me quedo con dicho patron y se lo asigno a x
                    x = transpose(xp(mu, :));

                    % Calculo la salida del perceptron ante la entrada x
                    y(mu) = signo(transpose(w) * x);

                    % Actaulizacion de w, aprendizaje.
                    w = w + constAprendizaje * (yd(mu) - y(mu)) * x;
                    error = computarError(y, yd);

                    cantPasos = cantPasos + 1;

                end

                pasos(i) = cantPasos;
            end

            iteraciones(iterConst) = mean(pasos);
            iterConst = iterConst + 1;

        end

        hold on;
        plot(constantes, iteraciones, 'linewidth', 3);
    end
    xlabel("Constante de aprendizaje"); ylabel("Cantidad de pasos");
    legend(legendTipoDeInit);
    set(gca,'fontsize', 12);
    set(grafico2,'PaperSize',[20 10]); %set the paper size to what you want  
    %print(grafico2,'Resultados/Ejercicio1/constAprVSiters','-dpdf') % then print it

    %scatter(constantes, iteraciones, 100, 'filled', 'black')
end

%% Perceptron simple para funcion AND de 4 entradas 
% Todas las salidas son -1 excepto para 1, 1, 1, 1

n = 5;
p = 2 ^ (n - 1);

xp = ones(p, n); % Patrones 
s = '0000';
for i = 1:p
    xp(i, 2:n) = s - '0';
    display(s);
    s = dec2bin(bin2dec(s) + bin2dec('0001'), 4);
end

yd = -1 * ones(p, 1); yd(p) = 1; % Salida esperada

% Inicio w al azar. Como son 3 entradas w es de tamaño (3,1). w =  w0,w1, w2
%w = rand(n, 1);
w = zeros(n, 1);
error = 1;
y = zeros(p, 1);
cantIter = 0;
constAprendizaje = 0.95;

while (error > 0)
    
    % Elijo un patron al azar
    mu = randi(p);
    disp("Patron " + mu);
    
    % Me quedo con dicho patron y se lo asigno a x
    x = transpose(xp(mu, :));
    disp("Patron actual: " + x);
    disp("w actual: " + w);
    
    % Calculo la salida del perceptron ante la entrada x
    y(mu) = signo(transpose(w) * x);
    
    % Actaulizacion de w, aprendizaje.
    w = w + constAprendizaje * (yd(mu) - y(mu)) * x;
    
    error = computarError(y, yd);
    disp("Error " + error);
    
    cantIter = cantIter + 1;
    
end

w0 = w(1); w1 = w(2); w2 = w(3); w3 = w(4); w4 = w(5);

display("Iteraciones: " + cantIter);
display(w);
display("y: " + y); display("yd: " + yd);

%% Perceptron Simple para una OR de 2 entradas

n = 3; % Cantidad de entradas
p = 4; % Canntidad de patrones


% Patrones. tamaño (4, 3) 4 patrones, 3 entradas
xp = [1, -1, -1; 1, -1, 1; 1, 1, -1; 1, 1, 1]; % Entrada 0 es 1
%xp = transpose(xp);
% Salida deseada (clasificacion 1 o -1)
yd = [-1; 1; 1; 1];

% Inicio w al azar. Como son 3 entradas w es de tamaño (3,1). w =  w0,w1, w2
%w = rand(n, 1);
w = zeros(n, 1);
error = 1;
y = zeros(p, 1);
cantIter = 0;
constAprendizaje = 0.9;

while (error > 0)
    
    % Elijo un patron al azar
    mu = randi(p);
    disp("Patron " + mu);
    
    % Me quedo con dicho patron y se lo asigno a x
    x = transpose(xp(mu, :));
    disp("Patron actual: " + x);
    disp("w actual: " + w);
    
    % Calculo la salida del perceptron ante la entrada x
    y(mu) = signo(transpose(w) * x);
    
    % Actaulizacion de w, aprendizaje.
    w = w + constAprendizaje * (yd(mu) - y(mu)) * x;
    
    error = computarError(y, yd);
    disp("Error " + error);
    
    cantIter = cantIter + 1;
    
end

w0 = w(1); w1 = w(2); w2 = w(3);
x1 = -2:0.1:2;
x2 = -w0 / w2 - x1 .* (w1 / w2);

display("Iteraciones: " + cantIter);
display(w);
display("y: " + y); display("yd: " + yd);

graficoOr = figure();
plot(x1, x2, 'linewidth', 3, 'color', 'red')
%fplot(@(t) -w0/w2 - t*w1/w2);
hold on;

scatter(xp(:, 2), xp(:, 3), 80, 'filled', 'black')
xlim([-2, 2]); ylim([-2, 2]); 
gridxy(0,0, 'Linestyle',':'); 
gridxy(1,1, 'Linestyle',':'); 
gridxy(-1,-1,'Linestyle',':');
title("Perceptron Simple OR 2 entradas");
set(gca,'fontsize', 14);
set(graficoOr,'PaperSize',[20 10]); %set the paper size to what you want  
%print(graficoOr,'Resultados/Ejercicio1/or2entradas','-dpdf') % then print it
%% Perceptron simple para funcion OR de 4 entradas 
% Todas las salidas son -1 excepto para 1, 1, 1, 1

n = 5;
p = 2 ^ (n - 1);

xp = ones(p, n); % Patrones 
s = '0000';
for i = 1:p
    xp(i, 2:n) = s - '0';
    display(s);
    s = dec2bin(bin2dec(s) + bin2dec('0001'), 4);
end

yd = ones(p, 1); yd(1) = -1; % Salida esperada

% Inicio w al azar. Como son 3 entradas w es de tamaño (3,1). w =  w0,w1, w2
%w = rand(n, 1);
w = zeros(n, 1);
error = 1;
y = zeros(p, 1);
cantIter = 0;
constAprendizaje = 0.95;

while (error > 0)
    
    % Elijo un patron al azar
    mu = randi(p);
    disp("Patron " + mu);
    
    % Me quedo con dicho patron y se lo asigno a x
    x = transpose(xp(mu, :));
    disp("Patron actual: " + x);
    disp("w actual: " + w);
    
    % Calculo la salida del perceptron ante la entrada x
    y(mu) = signo(transpose(w) * x);
    
    % Actaulizacion de w, aprendizaje.
    w = w + constAprendizaje * (yd(mu) - y(mu)) * x;
    
    error = computarError(y, yd);
    disp("Error " + error);
    
    cantIter = cantIter + 1;
    
end

w0 = w(1); w1 = w(2); w2 = w(3); w3 = w(4); w4 = w(5);

display("Iteraciones: " + cantIter);
display(w);
display("y: " + y); display("yd: " + yd);


