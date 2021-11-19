close all; clear all; clc;
rng('shuffle');

%% Modelo de Ising en 1D

N = 100; % Cantidad de dipolos 
temperaturas = [0.1:0.05:3];

% Genero la matriz W de conexiones entre dipolos
% w_ij = 1 si i y j son adyacentes, si no, vale 0
W = zeros(N, N);
for i = 1:N
    for j = 1:N
        if i == j
            W(i, j) = 0; % Conexion con si mismo es 0
        end
        if i > 1 && i < N 
            if j == i + 1 || j == i - 1
                W(i, j) = 1; % Conexion adyacente
            end
         % Condiciones de borde
        elseif i == 1 % Primer dipolo
            W(i, i + 1) = 1;
        else % Ultimo dipolo
            W(i, i - 1) = 1;
        end
    end
end
W = zeros(N) + diag(ones(N-1, 1),1) + diag(ones(N-1,1),-1);
% Inicialmente todos los dipolos estan magnetizados hacia arria (+1)
% Pruebo con una temperatura inicial T(1)
% - Defino recorrido = randperm(100); , defino vector de energias y el S inicial

% 1 - para una temperatura dada 
% 2 - calculo la energia del sistema
% 3 - elijo un spin de la lista recorrido y lo cambio
% 4 - calculo el deltaH
% 5 - acepto o no el cambio y guardo el nuevo <S>
% 6 - repito para 2-5 para todos los spines
% 7 - repito al menos 5 o 10 veces 2-6
% 8 - calculo la media <S>
% 9 - subo la temperatura y repito 2-8
% grafico <S> en func de la tempe

cantIter = 20; % Cantidad de iteraciones para metropolis.
S = ones(N, 1); % Arranco con todos los dipolos en +1
%S = 2 .* randi([0, 1], N, 1) - 1;
means_S = zeros(1, cantIter);
mean_S = zeros(1, length(temperaturas));
var_S = zeros(1, length(temperaturas)); 
hExt = 0; 
k = 1; 
idx = 1;
prevH = 0; 
afterH = 0;
deltaH = 0;

for T = temperaturas % Aumento la temp.

    for sweep = 1:cantIter % Repito el algoritmo Metropolis para sacar una media y varianza
        recorrido = randperm(100); % Lista al azar de indices de S para ir cambiando el estado de cada dipolo.

        % Metropolis
        for i = 1:N % iterio para todos los spines al azar
            
            % Calculo la energia del sistema 
            prevH = -1/2 * S' * W  * S - hExt * sum(S);
            
            % Cambio el estado de un spin
            S(recorrido(i)) = - S(recorrido(i));
            
            % Calculo la nueva energia
            afterH = -1/2 * S' * W  * S - hExt * sum(S);
            
            % Cambio de energia
            deltaH = afterH - prevH;

            % Acepto el cambio con probabilidad, si el cambio de energia
            % fue mayor a 0.
            if deltaH > 0
                % Genero nro de 0 a 1
                r = rand;
                prob = exp( - deltaH / (k * T) );
                
                if r > prob
                    % No acepto el cambio, por ende vuelvo a cambiarlo
                    S(recorrido(i)) = - S(recorrido(i));
                end
            end
            
        end
        
        % Calculko la media de S
        means_S(sweep) = sum(S) / N;
    end
    
    mean_S(idx) = sum(means_S) / cantIter;
    var_S(idx) = sqrt(var(means_S));
    idx = idx + 1;
end

h1 = figure()
%errorbar( temperaturas, mean_S, var_S)
scatter(temperaturas, mean_S, 'filled', 'black')
grid('on')
title("Ising en 1D");
legend("<S>");
grid('on');
xlabel('Temperatura');
ylabel('<S>');
set(gca,'fontsize', 10);
set(gca, 'XDir','reverse')
set(h1,'PaperSize',[20 10]); %set the paper size to what you want  
print(h1,'Resultados/Ejercicio4/1d','-dpdf') % then print it
%ylim([-0.2, 1]);

%% Modelo de Ising para 2D
n = 20;
N = n*n; % Cantidad de dipolos 
W = generarConexiones(N, n);
temperaturas = [0.1:0.05:4];
%temperaturas = generarTemp(80, 0.9);
cantIter = 15; % Cantidad de iteraciones para metropolis.
S = ones(N, 1); % Arranco con todos los dipolos en +1
%S = 2 .* randi([0, 1], N, 1) - 1;
means_S = zeros(1, cantIter);
mean_S = zeros(1, length(temperaturas));
var_S = zeros(1, length(temperaturas)); 
hExt = 0; 
k = 1; 
idx = 1;
prevH = 0; 
afterH = 0;
deltaH = 0;

for T = temperaturas 

    for iter = 1:cantIter % repito Metropolis cantIter veces
        recorrido = randperm(N); % Lista al azar de los indices de S

        for i = 1:N % Itero sobre todos los dipolos 
            
            % Calculo la energia del sistema 
            prevH = -1/2 * S' * W  * S - hExt * sum(S);
            
            % Cambio el estado de un spin
            S(recorrido(i)) = - S(recorrido(i));
            
            % Calculo la nueva energia
            afterH = -1/2 * S' * W  * S - hExt * sum(S);
            
            % Cambio de energia
            deltaH = afterH - prevH;

            % Acepto el cambio con probabilidad, si el cambio de energia
            % fue mayor a 0.
            if deltaH > 0
                % Genero nro de 0 a 1
                r = rand;
                prob = exp( - deltaH / (k * T) );
                
                if r > prob
                    % No acepto el cambio, por ende vuelvo al valor
                    % original
                    S(recorrido(i)) = - S(recorrido(i));
                end
            end

        end
        
        means_S(iter) = sum(S) / N;
        
    end
    
    mean_S(idx) = sum(means_S) / cantIter;
    var_S(idx) = sqrt(var(means_S));
    
    idx = idx + 1;
end
h2=figure()
%errorbar( temperaturas, mean_S, var_S)
scatter(temperaturas, mean_S, 'filled', 'black')
grid('on')
title("Ising en 2D");
legend("<S>");
grid('on');
xlabel('Temperatura');
ylabel('<S>');
set(gca,'fontsize', 10);
%ylim([-0.2, 1]);
set(gca, 'XDir','reverse')
%set(h2,'PaperSize',[20 10]); %set the paper size to what you want  
%print(h2,'Resultados/Ejercicio4/2d','-dpdf') % then print it

%hold('on')
%       a =   2.912e+06  ;
%       b =      -21.63 ;
%       c =   2.873e+06 ;
%fittingCurve = a./(exp(-b.*temperaturas) + c);

%plot(temperaturas, fittingCurve, 'linewidth', 1.5, 'color', 'red');

%% Generar matriz de dipolos en 2D
% N n[umero de dipolos
% n dimension de la matriz
% por ejemplo 20*20 dipolos = N y n = 20
function Y = generarConexiones(N,n)
    aux1 = ones(N - 1, 1);
    for i = 1:(N - 1)
       if mod(i, n) == 0
           aux1(i) = 0;
       end
    end    
    aux2 = ones(N - n, 1);
    Y = diag(aux1, 1) + diag(aux1, -1) + diag(aux2, n) + diag(aux2, -n);
end
