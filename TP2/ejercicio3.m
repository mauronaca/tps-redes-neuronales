
%% a) Perceptron multicapa que aprende la funcion XOR de 2 entradas
% La red tendrï¿½ 3 capas: La capa de entrada x, la capa oculta que tendrï¿½ 2
% neuronas dado que la dimensiï¿½n de entrada es de 2, tiene sentido que
% esta capa tenga mï¿½nimo 2 neuronas. La tercera capa de salida serï¿½ una
% neurona ya que la salida es 1  o -1
%
clear all; clc; close all;
rng('shuffle');


N_0 = 2 + 1; % 2 entradas de la funciï¿½n + 1 de bias
N_1 = 3; % #neuronas de la capa oculta 
N_2 = 1; % #neuronas de la capa de salida
N = {N_0, N_1, N_2}; % #neuronas de cada capa.

p = 4; % #patrones
n = 2; % #capas

% patrones de entrada, cada fila es un patron, cada columna una entrada.
xp = [1, -1, -1; 1, -1, 1; 1, 1, -1; 1, 1, 1];

yd = [-1; 1; 1; -1]; % Salida deseada para cada patron -> tamaï¿½o de (p, 1).
Vd = yd;
y = zeros(p, 1); % Defino el vector salida del perceptron.

h = {zeros(N_1, 1), zeros(N_2, 1)}; % Estados de las neuronas para cada capa
V = {zeros(N_0, 1), tanh(h{1}), tanh(h{2})}; % Salida de cada neurona para cada capa. La primer
% componente es la entrada y a partir del segundo es la capa 1, y luego la
% capa 2.
delta = {zeros(N_1, 1), zeros(N_2, 1)}; % Errores de c/neuronas en cada capa.

constAprendizaje = 0.9;
error = 1;


% Inicializo los pesos de la red:
% Para cada matriz de pesos, w_ij donde i es la neurona de la capa actual,
% y j de la capa anterior.
w = {randi([-10, 10], N_1, N_0)/10, randi([-10,10], N_2, N_1)/10}; % Pesos para cada capa
bias = {rand(N_1, 1), rand(N_2, 1)};

%delta_w = zeros(1,p);
delta_w = {zeros(N_1, N_0), zeros(N_2, N_1)};

% iteraciones
iter = 0; iter_max = 1000;
tol = 0.01;

% Vector de errores
errores = zeros(1, iter_max);

plot_xor2 = figure();
while (error > tol) && (iter < iter_max)
    delta_w = {zeros(N_1, N_0), zeros(N_2, N_1)};
    
    for mu = randperm(p)
        % Itero sobre cada patrï¿½n mu.
        x = (xp(mu, :))';
        %disp("Patron "+mu+":"); %disp(x);
        V{1} = x;
        
        %
        % Obtengo la salida real:
        % Itero por cada capa, y en cada una de ellas itero sobre cada
        % neurona. Entonces calculo primero el estado 'hi' de la neurona i
        % haciendo la sumatoria de los w_ij * Vj, donde Vj es la salida de
        % la neurona de la capa anterior. (para la capa 1, Vj serï¿½ la entrada x)
        for m = 1:n
            h{m} = w{m} * V{m} + bias{m};
            V{m+1} = tanh(h{m});
        end
        
        %
        % Computa el delta error empezando por la capa de salida hasta la
        % primera capa
        %
        % Primer el error de la capa de salida, entre la salida y la
        % deseada.
        delta{n} = (1 - tanh(h{n})^2) * (Vd(mu) - V{n+1});
        % Ahora se propaga el error hacia atras.
        for m=n-1:-1:1
            % Itero por cada neurona:
            %for i=1:size(h{m}, 1)
            %    delta{m}(i) = (1 - tanh(h{m}(i)^2)) * (transpose(w{m+1}(:, i)) * delta{m+1});
            %end
            delta{m} = ( 1 - tanh(h{m}).^2 ) .* ( w{m+1} * delta{m+1} )';
        end
        
        %
        % Computo delta_w para cada patron:
        % Voy a ir acumulando los resultados para cada delta_w_ij en vez de
        % guardar los delta_w de cada patron
        for m=1:n 
            delta_w{m} = delta_w{m} + constAprendizaje * delta{m} * V{m}';
        end
           
    end
    
    % Acvtualizo los w utilizando los delta_w
    for m = 1:n
        w{m} = w{m} + delta_w{m};
    end
    
    % Computo la salida para cada patron
    for mu = 1:p
        y(mu) = computarSalida(transpose(xp(mu, :)), w, bias);
    end
    
    % Computo el error
    error = computarError(y, yd);
    disp("Iteracion " + iter + " y error " + error);
    iter = iter + 1;
    errores(iter) = error;

end

if (iter == iter_max)
    disp("No convergió, error " + error);
else
    disp("Cantidad de iteraciones " + iter + " y error " + error);
    for mu = 1:p
        x = transpose(xp(mu,:));
        %disp(x);
        salida = computarSalida(x, w, bias);
        %disp(salida);
        plot(1:iter, errores(1,1:iter), 'linewidth', 3, 'color', 'black');
        xlim([1,iter]);
        grid('on');
        set(gca,'fontsize', 14);
        ylabel('Error'); xlabel('Iteración')
        set(plot_xor2,'PaperSize',[20 10]); %set the paper size to what you want  
        print(plot_xor2,'Resultados/Ejercicio3/xor2','-dpdf') % then print it
    end
end

%%% Punto C:
% Cambiar pesos y observar mibimos y mesetas!
%
cant_cambios = 1000;
peso = zeros(1, cant_cambios);
errores_cambios = zeros(1, 2*cant_cambios);
errores_mu = {zeros(1, 2*cant_cambios),zeros(1, 2*cant_cambios),zeros(1, 2*cant_cambios),zeros(1, 2*cant_cambios)};
iters = 1;

% Elijo un peso al azarrr
m1 = randi([1,2]);
i1 = randi([1,size(w{m1},1)]);
j1 = randi([1,size(w{m1},2)]);
wOriginal1 = w{m1}(i1, j1); 

m2 = randi([1,2]);
i2 = randi([1,size(w{m2},1)]);
j2 = randi([1,size(w{m2},2)]);
wOriginal2 = w{m2}(i2, j2);

wNuevos = (-cant_cambios:1:cant_cambios) / 100;

for wNuevo = wNuevos
    % Cambio el peso al azar   
    w{m1}(i1, j1) = wNuevo;  
    %w{m2}(i2, j2) = wNuevo;  
    
    for mu = 1:p
        x = xp(mu, :)';
        y(mu) = computarSalida(x, w, bias);
        errores_mu{mu}(iters) = computarError(y(mu), yd(mu));
    end
    
    errores_cambios(iters) = computarError(y, yd);
    iters = iters + 1;
end
plot_variar_pesos = figure();
plot(wNuevos, errores_cambios, 'linewidth', 3, 'color', 'black');
hold on;
scatter(wOriginal1, error, 140 , 'filled' , 'red')
%scatter(wOriginal2, error, 140 , 'filled' , 'blue')
legend('Error', strcat('Peso utilizado w',num2str(i1),num2str(j1)));
xlabel(strcat('Variación del peso w',num2str(i1),num2str(j1)));
ylabel('Error');
set(gca,'fontsize', 12);
grid('on');
set(plot_variar_pesos, 'PaperSize',[20 10]); %set the paper size to what you want  
print(plot_variar_pesos,'Resultados/Ejercicio3/puntoC','-dpdf') % then print it

plot_variar_pesos_patron = figure();
colors = ["red" "magenta" "blue" "black"];
for plots = 1:p
    hold on;
    plot(wNuevos, errores_mu{plots}, 'linewidth', 1.5, 'color',colors(plots));
end
hold on;
scatter(wOriginal1, error, 140 , 'filled' , 'red')
legend('Error y1','Error y2','Error y3','Error y4', strcat('Peso utilizado w',num2str(i1),num2str(j1)));
xlabel(strcat('Variación del peso $w_$',num2str(i1),num2str(j1)));
ylabel('Error');
set(gca,'fontsize', 12);
grid('on');
set(plot_variar_pesos_patron, 'PaperSize',[20 10]); %set the paper size to what you want  
print(plot_variar_pesos_patron,'Resultados/Ejercicio3/puntoD','-dpdf') % then print it

pause;
%% Red multicapa que aprende la funcion XOR de 4 entradas

% La red tendrá 3 capas, la entrada, la capa oculta y la capa de salida. La
% entrada será de dimensión 4 + 1 de bias. La capa oculta es logico pensar
% que tenga misma cantidad de neuronas que la entrada, así que tendrá 5
% neuronas. La capa de salida 1 neurona ya que la salida será 1 o -1.
% 
clear all; clc; close all;
rng('shuffle');

n = 4 + 1; % #entradas + bias
p = 2^4;
N = [5, 1]; % #neuronas por capa
M = 2; % #capas

% Patrones; 
xp = ones(p, n); % Cada fila es un patron distinto. Y cada columna corresponde a una entrada.
s = '0000';
for i = 1:p
    xp(i, 2:n) = s - '0';
    s = dec2bin(bin2dec(s) + bin2dec('0001'), 4);
end
xp(:, 1) = 1; % Bias!
xp = 2 * xp - 1;

yd = ones(p, 1); yd(1) = -1; yd(p) = -1; % Salida deseada
Vd = yd;
y = zeros(p, 1);

h = {zeros(N(1), 1), zeros(N(2), 1)}; % Estados de las neuronas en c/capa
V = {zeros(n, 1), tanh(h{1}), tanh(h{2})}; % Salida de c/neurona para c/capa. El primer elemento es la entrada y la
% inicializo con el primer patron.

sigma = {zeros(N(1), 1), zeros(N(2), 1)}; % Errores de c/neuronas en cada capa. Es sigma, no delta

% Inicializo los pesos:
delta_w = {zeros(N(1),n), zeros(N(2), N(1))}; % 
%w = {rand(N(1), n) , -rand(N(2), N(1)) };
w = {randi([-10,10], N(1), n) / 10 , randi([-10,10], N(2), N(1))/10};
bias = {rand(N(1), 1), rand(N(2), 1)};

error = 1;
constAprendizaje = 0.1;
iterMax = 1000;
iter = 0;
tol = 0.1;

% Vector de errores
errores = zeros(1, iterMax);

plot_xor4 = figure();
while(error > tol) && (iter < iterMax)
    
    delta_w = {zeros(N(1),n), zeros(N(2), N(1))};
    
    % Itero para cada patron de manera aleatoria:
    for mu = 1:p
        %disp("Patron: " + mu);
        x = xp(mu, :)';
        V{1} = x; % La entrada de la red.
        
        % Obtengo la salida real, a partir de presentar el patron de
        % entrada a la red:
        for m = 1:M
            h{m} = w{m} * V{m} + bias{m};
            V{m+1} = tanh(h{m});
        end
        
        % Backpropagation: Calculo el error desde la capa de salida hacia
        % la de entrada:
        % Primero se calcula el error entre la salida y la entrada en la
        % ultima capa:
        sigma{M} = (1 - tanh(h{M})^2) * (Vd(mu) - V{M+1});
        % las demas capas:
        for m  = M-1:-1:1
            %disp(sigma{m+1})
            %disp(w{m+1}')
            %disp(h{m})
            
            sigma{m} = ( 1 - tanh(h{m}).^2 ) .* ( w{m+1} * sigma{m+1} )';
            %itero para cada neurona
            %for i = randperm(N(m))
         
             %   sigma{m}(i) = ( 1 - tanh(h{m}(i))^2 ) * ( w{m+1}(:, i)' * sigma{m+1} );
            %end
        end
    
        % Calculo los delta w
        for m=1:M
            %for i = randperm(size(w{m},1))
            %    for j = randperm(size(w{m},2))
            %        delta_w{m}(i, j) = delta_w{m}(i, j) + constAprendizaje * sigma{m}(i) * V{m}(j); 
            %    end
            %end
            delta_w{m} = delta_w{m} + constAprendizaje * sigma{m} * V{m}';
        end
  
        
    end
    
    % Actualizo los w
    for m = 1:M
        w{m} = w{m} + delta_w{m};
     
    end
    %disp("Pesos:");disp(w{1});disp(w{2});
    % Computo la salida
    for mu = 1:p
        y(mu) = computarSalida(xp(mu, :)', w, bias);
        %error = error + 0.5*(yd(mu)-y(mu))^2;
    end
    
    % Computo el error:
    error = computarError(y, yd);
    disp("Iteracion " + iter + " y error " + error);
    iter = iter + 1;
    errores(iter) = error;
end

if (iter == iterMax)
    disp("No convergió, error " + error);
else
    disp("Cantidad de iteraciones " + iter + " y error " + error);
    for mu = 1:p
        x = transpose(xp(mu,:));
        %disp(x);
        salida = computarSalida(x, w, bias);
        %disp(salida);
        plot(1:iter, errores(1,1:iter), 'linewidth', 3, 'color', 'black');
        xlim([1,iter]);
        grid('on');
        set(gca,'fontsize', 14);
        ylabel('Error'); xlabel('Iteración')
        set(plot_xor4,'PaperSize',[20 10]); %set the paper size to what you want  
        %print(plot_xor4,'Resultados/Ejercicio3/xor4','-dpdf') % then print it
    end
end
