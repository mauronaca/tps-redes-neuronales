
%% Red con aprendizaje genetico
% Se busca contar cuantas generaciones se requieren antes de que converja,
% para distintas variaciones de eta, crossover, y tamaño de poblacion
clc; clear all; close all;

sigma = rand() * 0.2; % tasa de mutacion
p_cross = 0.35; % proba. de crossover
Nind = 50; % #individuos por generacion

n = 2 + 1; % #entradas
p = 4; % #patrones
Nneu = 3; % #neuronas capa oculta
Nsal = 1; % #neuronas capa salida

% patrones
x0 = ones(p,1);
x1 = [-1; -1; 1; 1];
x2 = [-1; 1; -1; 1];
X = [x0 x1 x2];

yd = [-1; 1; 1; -1]; % salida deseada


% Inicializo los pesos de cada individuo. Cada celda es un individuo...
w1_poblacion = cell(1, Nind);
w2_poblacion = cell(1, Nind);

b1 = rand(Nneu, 1) - 0.5;
b2 = rand() - 0.5;

for i = 1:Nind
    w1_poblacion{i} = rand(Nneu, n) - 0.5;
    w2_poblacion{i} = rand(Nsal, Nneu) - 0.5;
end


% se calcula el fitness inicial
vec_fitness = zeros(Nind , 1); 
error = 0;
tolerancia = 0;
fitness_d = 1 - tolerancia /4; % El fitness deseado, queremos maximizarlo
fitness_elite = 0;
fitness_max = fitness_elite;

for individuo = 1: Nind
    error = 0;
    
    for mu = 1:p
        x = X(mu, :)';
        
        h1 = w1_poblacion{individuo} * x + b1;
        V1 = signo(h1);
        
        h2 = w2_poblacion{individuo} * V1 + b2;
        V2 = signo(h2);
        
        error = error + (yd(mu) - V2)^2;
    end
    
    error = error / p;
    vec_fitness(individuo) = 1 - error / 4;
end

[fitness_elite, elite_idx] = max(vec_fitness);

if fitness_elite > vec_fitness(end)
    fitness_max = [fitness_max fitness_elite];
else
    fitness_elite = fitness_max(end);
    fitness_max = [fitness_max fitness_elite];
end

% aprendizaje de la red:
iter_max = 1000;
iter = 0;

while (fitness_elite < fitness_d) && (iter < iter_max)
    
    % me guardo los indiviudos de elite
    w1_elite = w1_poblacion{elite_idx};
    w2_elite = w2_poblacion{elite_idx};
    
    p_repr = vec_fitness .* 1 / sum(vec_fitness); % proba de reproduccion
    
    idx_repr = randsample(1:Nind, Nind, true, p_repr);
    
    % inicializo la poblacion de reproducidos
    w1_poblacion_repr = cell(1, Nind);
    w2_poblacion_repr = cell(1, Nind);
    for i = 1:Nind
        w1_poblacion_repr{i} = zeros(Nneu, n) - 0.5;
        w2_poblacion_repr{i} = zeros(Nsal, Nneu) - 0.5;
    end
    
    for individuo1 = 1:(Nind - 1)
        for individuo2 = idx_repr
            w1_poblacion_repr{individuo1} = w1_poblacion{individuo2};
            w2_poblacion_repr{individuo1} = w2_poblacion{individuo2};
            
            % crossover
            r = rand();
            if r < p_cross
                % se acepta el cruce entre individuos
                % con uno distinto al individuo actual
                individuo3 = idx_repr(randi(Nind));
                while(individuo3 == individuo2)
                    individuo3 = idx_repr(randi(Nind));
                end
                % Se hace la cruza:
                w1_poblacion_repr{individuo1} = w1_poblacion{individuo3};
                w2_poblacion_repr{individuo1} = w2_poblacion{individuo3};
            end
            
        end
    end
    
    % mutacion
    w1_delta = cell(1, Nind);
    w2_delta = cell(1, Nind);
    for i=1:Nind
        w1_delta{i} = generar_w_mutacion(Nneu, n, sigma);
        w2_delta{i} = generar_w_mutacion(Nsal, Nneu, sigma);
    end    
    
    for i=1:Nind
        w1_poblacion_repr{i} = w1_poblacion_repr{i} + w1_delta{i};
        w2_poblacion_repr{i} = w2_poblacion_repr{i} + w2_delta{i};
    end
    w1_poblacion_repr{Nind} = w1_elite;
    w2_poblacion_repr{Nind} = w2_elite;
    
    w1_poblacion = w1_poblacion_repr;
    w2_poblacion = w2_poblacion_repr;
    
    
    % recalcula el fitness 
    vec_fitness = zeros(Nind , 1);
    error = 0;
    
    for individuo = 1: Nind
        error = 0;
        
        for mu = 1:p
            x = X(mu, :)';

            h1 = w1_poblacion{individuo} * x + b1;
            V1 = signo(h1);

            h2 = w2_poblacion{individuo} * V1 + b2;
            V2 = signo(h2);

            error = error + (yd(mu) - V2)^2;
        end

        error = error / p;
        vec_fitness(individuo) = 1 - error / 4;
    end

    [fitness_elite, elite_idx] = max(vec_fitness);

    if fitness_elite > vec_fitness(end)
        fitness_max = [fitness_max fitness_elite];
    else
        fitness_elite = fitness_max(end);
        fitness_max = [fitness_max fitness_elite];
    end
    
    iter = iter + 1;
    display("Iteracion: " + iter);
end

plot_a = figure();
plot(1:size(fitness_max, 2), fitness_max, 'linewidth', 3, 'color', 'black');
xlabel('Generación'); ylabel('Fitness'); 
legend("sigma = " + sigma + newline + "prob. crossover = " + p_cross + newline + "cant. individuos = " + Nind + newline + "iteraciones = " + iter, 'location', 'northwest');
set(gca,'fontsize', 12); xlim([1,size(fitness_max, 2)]);
set(plot_a, 'PaperSize', [20 10]); 
print(plot_a,'Resultados/Ejercicio5/fitnessVSgen.pdf','-dpdf');

%% Generacion vs Sigma 
% Ahora se hacen varias realizaciones para ver la cantidad de generaciones que le lleva converger en funcion 
% de la tasa de ,mutacion
clc; clear all; close all;

Nind = 10; % #individuos
p_cross = 0.35; % proba. de crossover

n = 2 + 1; % #entradas
p = 4; % #patrones
Nneu = 3; % #neuronas capa oculta
Nsal = 1; % #neuronas capa salida

% patrones
x0 = ones(p,1);
x1 = [-1; -1; 1; 1];
x2 = [-1; 1; -1; 1];
X = [x0 x1 x2];

yd = [-1; 1; 1; -1]; % salida deseada

vec_sigma = 0.05:.05:0.95; %tasas de mutaciones para realizaciones distitas
cant_sigmas = size(vec_sigma, 2);

% Inicializo los pesos de cada individuo. Guardo la inicializacion para
% cada iteracion
w1_poblacion_inicial = cell(1, Nind); % capa oculta
w2_poblacion_inicial = cell(1, Nind); % capa salida
b1_inicial = rand(Nneu, 1) - 0.5; % biasess
b2_inicial = rand() - 0.5;

for i = 1:Nind
    w1_poblacion_inicial{i} = rand(Nneu, n) - 0.5;
    w2_poblacion_inicial{i} = rand(Nsal, Nneu) - 0.5;
end

% iteraciones / Generaciones
gen_max = 500; % si nunca converge, corta por aca 
generaciones = zeros(size(vec_sigma));

cant_iteraciones = 100;


for sigma = vec_sigma
    % itero para distintas tasas de mutacion
    gen = 0;
    idx = round(sigma * cant_sigmas);
    display("Sigma: "+ sigma)
    for realizacion = 1:cant_iteraciones
        % entreno la red!
        
        % inicializo los pesos, siempre con el mismo valor random generado
        % al ppio.
        w1_poblacion = w1_poblacion_inicial;
        w2_poblacion = w2_poblacion_inicial;
        % bias inicializo
        b1 = b1_inicial;
        b2 = b2_inicial;
        
        error = 1;
        tolerancia = 0;
        fitness_d = 1 - tolerancia / 4; % fitness target
        fitness_elite = 0;
        fitness_max = fitness_elite;
        
        % calculo inicialmente el fitness:
        for individuo = 1: Nind
            error = 0;

            for mu = 1:p
                x = X(mu, :)';

                h1 = w1_poblacion{individuo} * x + b1;
                V1 = signo(h1);

                h2 = w2_poblacion{individuo} * V1 + b2;
                V2 = signo(h2);

                error = error + (yd(mu) - V2)^2;
            end

            error = error / p;
            vec_fitness(individuo) = 1 - error / 4;
        end

        [fitness_elite, elite_idx] = max(vec_fitness);

        if fitness_elite > vec_fitness(end)
            fitness_max = [fitness_max fitness_elite];
        else
            fitness_elite = fitness_max(end);
            fitness_max = [fitness_max fitness_elite];
        end
        
        % aprendizaje de la red:
        gen = 1;
        while (fitness_elite < fitness_d) && (gen < gen_max)

            % me guardo los indiviudos de elite
            w1_elite = w1_poblacion{elite_idx};
            w2_elite = w2_poblacion{elite_idx};

            p_repr = vec_fitness .* 1 / sum(vec_fitness); % proba de reproduccion

            idx_repr = randsample(1:Nind, Nind, true, p_repr);

            % inicializo la poblacion de reproducidos
            w1_poblacion_repr = cell(1, Nind);
            w2_poblacion_repr = cell(1, Nind);
            for i = 1:Nind
                w1_poblacion_repr{i} = zeros(Nneu, n) - 0.5;
                w2_poblacion_repr{i} = zeros(Nsal, Nneu) - 0.5;
            end

            for individuo1 = 1:(Nind - 1)
                for individuo2 = idx_repr
                    w1_poblacion_repr{individuo1} = w1_poblacion{individuo2};
                    w2_poblacion_repr{individuo1} = w2_poblacion{individuo2};

                    % crossover
                    r = rand();
                    if r < p_cross
                        % se acepta el cruce entre individuos
                        % con uno distinto al individuo actual
                        individuo3 = idx_repr(randi(Nind));
                        while(individuo3 == individuo2)
                            individuo3 = idx_repr(randi(Nind));
                        end
                        % Se hace la cruza:
                        w1_poblacion_repr{individuo1} = w1_poblacion{individuo3};
                        w2_poblacion_repr{individuo1} = w2_poblacion{individuo3};
                    end

                end
            end

            % mutacion
            w1_delta = cell(1, Nind);
            w2_delta = cell(1, Nind);
            for i=1:Nind
                w1_delta{i} = generar_w_mutacion(Nneu, n, sigma);
                w2_delta{i} = generar_w_mutacion(Nsal, Nneu, sigma);
            end    

            for i=1:Nind
                w1_poblacion_repr{i} = w1_poblacion_repr{i} + w1_delta{i};
                w2_poblacion_repr{i} = w2_poblacion_repr{i} + w2_delta{i};
            end
            w1_poblacion_repr{Nind} = w1_elite;
            w2_poblacion_repr{Nind} = w2_elite;

            w1_poblacion = w1_poblacion_repr;
            w2_poblacion = w2_poblacion_repr;


            % recalcula el fitness 
            vec_fitness = zeros(Nind , 1);
            error = 0;

            for individuo = 1: Nind
                error = 0;

                for mu = 1:p
                    x = X(mu, :)';

                    h1 = w1_poblacion{individuo} * x + b1;
                    V1 = signo(h1);

                    h2 = w2_poblacion{individuo} * V1 + b2;
                    V2 = signo(h2);

                    error = error + (yd(mu) - V2)^2;
                end

                error = error / p;
                vec_fitness(individuo) = 1 - error / 4;
            end

            [fitness_elite, elite_idx] = max(vec_fitness);

            if fitness_elite > vec_fitness(end)
                fitness_max = [fitness_max fitness_elite];
            else
                fitness_elite = fitness_max(end);
                fitness_max = [fitness_max fitness_elite];
            end

            gen = gen + 1;
            %display("Generacion: " + gen);
        end
        generaciones(idx) = generaciones(idx) + gen;
    end
    
    generaciones(idx) = generaciones(idx) /  cant_iteraciones;
    % siguiente sigmas
end

plot_b = figure();
plot(vec_sigma, generaciones, 'linewidth', 3, 'color', 'black');
xlabel('Tasa de mutación'); ylabel('Generaciones'); 
legend("prob. crossover = " + p_cross + newline + "cant. individuos = " + Nind + newline + "cant. realizaciones = " + cant_iteraciones, 'location', 'northwest');
set(gca,'fontsize', 12); 
set(plot_b, 'PaperSize', [20 10]); 
print(plot_b,'Resultados/Ejercicio5/ptoB.pdf','-dpdf');



%% Generacion vs proba de crossover 
% #generaciones que le lleva converger en funcion de la proba de crossover.
clc; clear all; close all;

Nind = 10; % #individuos

n = 2 + 1; % #entradas
p = 4; % #patrones
Nneu = 3; % #neuronas capa oculta
Nsal = 1; % #neuronas capa salida

% patrones
x0 = ones(p,1);
x1 = [-1; -1; 1; 1];
x2 = [-1; 1; -1; 1];
X = [x0 x1 x2];

yd = [-1; 1; 1; -1]; % salida deseada

vec_pcross = 0.05:0.05:0.95; % proba. de crossover
cant_pcross = size(vec_pcross, 2);

sigma = 0.45; %tasas de mutaciones para realizaciones distitas

% Inicializo los pesos de cada individuo. Guardo la inicializacion para
% cada iteracion
w1_poblacion_inicial = cell(1, Nind); % capa oculta
w2_poblacion_inicial = cell(1, Nind); % capa salida
b1_inicial = rand(Nneu, 1) - 0.5; % biasess
b2_inicial = rand() - 0.5;

for i = 1:Nind
    w1_poblacion_inicial{i} = rand(Nneu, n) - 0.5;
    w2_poblacion_inicial{i} = rand(Nsal, Nneu) - 0.5;
end

% iteraciones / Generaciones
gen_max = 500; % si nunca converge, corta por aca 
generaciones = zeros(size(vec_pcross));

cant_iteraciones = 100;


for p_cross = vec_pcross
    % itero para distintas tasas de mutacion
    gen = 0;
    idx = round(p_cross * cant_pcross);
    display("Proba de crossover: "+ p_cross)
    for realizacion = 1:cant_iteraciones
        % entreno la red!
        
        % inicializo los pesos, siempre con el mismo valor random generado
        % al ppio.
        w1_poblacion = w1_poblacion_inicial;
        w2_poblacion = w2_poblacion_inicial;
        % bias inicializo
        b1 = b1_inicial;
        b2 = b2_inicial;
        
        error = 1;
        tolerancia = 0;
        fitness_d = 1 - tolerancia / 4; % fitness target
        fitness_elite = 0;
        fitness_max = fitness_elite;
        
        % calculo inicialmente el fitness:
        for individuo = 1: Nind
            error = 0;

            for mu = 1:p
                x = X(mu, :)';

                h1 = w1_poblacion{individuo} * x + b1;
                V1 = signo(h1);

                h2 = w2_poblacion{individuo} * V1 + b2;
                V2 = signo(h2);

                error = error + (yd(mu) - V2)^2;
            end

            error = error / p;
            vec_fitness(individuo) = 1 - error / 4;
        end

        [fitness_elite, elite_idx] = max(vec_fitness);

        if fitness_elite > vec_fitness(end)
            fitness_max = [fitness_max fitness_elite];
        else
            fitness_elite = fitness_max(end);
            fitness_max = [fitness_max fitness_elite];
        end
        
        % aprendizaje de la red:
        gen = 1;
        while (fitness_elite < fitness_d) && (gen < gen_max)

            % me guardo los indiviudos de elite
            w1_elite = w1_poblacion{elite_idx};
            w2_elite = w2_poblacion{elite_idx};

            p_repr = vec_fitness .* 1 / sum(vec_fitness); % proba de reproduccion

            idx_repr = randsample(1:Nind, Nind, true, p_repr);

            % inicializo la poblacion de reproducidos
            w1_poblacion_repr = cell(1, Nind);
            w2_poblacion_repr = cell(1, Nind);
            for i = 1:Nind
                w1_poblacion_repr{i} = zeros(Nneu, n) - 0.5;
                w2_poblacion_repr{i} = zeros(Nsal, Nneu) - 0.5;
            end

            for individuo1 = 1:(Nind - 1)
                for individuo2 = idx_repr
                    w1_poblacion_repr{individuo1} = w1_poblacion{individuo2};
                    w2_poblacion_repr{individuo1} = w2_poblacion{individuo2};

                    % crossover
                    r = rand();
                    if r < p_cross
                        % se acepta el cruce entre individuos
                        % con uno distinto al individuo actual
                        individuo3 = idx_repr(randi(Nind));
                        while(individuo3 == individuo2)
                            individuo3 = idx_repr(randi(Nind));
                        end
                        % Se hace la cruza:
                        w1_poblacion_repr{individuo1} = w1_poblacion{individuo3};
                        w2_poblacion_repr{individuo1} = w2_poblacion{individuo3};
                    end

                end
            end

            % mutacion
            w1_delta = cell(1, Nind);
            w2_delta = cell(1, Nind);
            for i=1:Nind
                w1_delta{i} = generar_w_mutacion(Nneu, n, sigma);
                w2_delta{i} = generar_w_mutacion(Nsal, Nneu, sigma);
            end    

            for i=1:Nind
                w1_poblacion_repr{i} = w1_poblacion_repr{i} + w1_delta{i};
                w2_poblacion_repr{i} = w2_poblacion_repr{i} + w2_delta{i};
            end
            w1_poblacion_repr{Nind} = w1_elite;
            w2_poblacion_repr{Nind} = w2_elite;

            w1_poblacion = w1_poblacion_repr;
            w2_poblacion = w2_poblacion_repr;


            % recalcula el fitness 
            vec_fitness = zeros(Nind , 1);
            error = 0;

            for individuo = 1: Nind
                error = 0;

                for mu = 1:p
                    x = X(mu, :)';

                    h1 = w1_poblacion{individuo} * x + b1;
                    V1 = signo(h1);

                    h2 = w2_poblacion{individuo} * V1 + b2;
                    V2 = signo(h2);

                    error = error + (yd(mu) - V2)^2;
                end

                error = error / p;
                vec_fitness(individuo) = 1 - error / 4;
            end

            [fitness_elite, elite_idx] = max(vec_fitness);

            if fitness_elite > vec_fitness(end)
                fitness_max = [fitness_max fitness_elite];
            else
                fitness_elite = fitness_max(end);
                fitness_max = [fitness_max fitness_elite];
            end

            gen = gen + 1;
            %display("Generacion: " + gen);
        end
        generaciones(idx) = generaciones(idx) + gen;
    end
    
    generaciones(idx) = generaciones(idx) /  cant_iteraciones;
    % siguiente sigmas
end


plot_c = figure();
plot(vec_pcross, generaciones, 'linewidth', 3, 'color', 'black');
xlabel('Probabilidad de Crossover'); ylabel('Generaciones'); 
legend("sigma = " + sigma + newline  + "cant. individuos = " + Nind , 'location', 'northwest');
set(gca,'fontsize', 12); 
set(plot_c, 'PaperSize', [20 10]); 
print(plot_c,'Resultados/Ejercicio5/genVSpcross.pdf','-dpdf');


%% Generaciones vs Cant. Individuos

clc; clear all; close all;


vec_individuos = 10:10:100; % #individuos
cant_ind = size(vec_individuos, 2);

n = 2 + 1; % #entradas
p = 4; % #patrones
Nneu = 3; % #neuronas capa oculta
Nsal = 1; % #neuronas capa salida

% patrones
x0 = ones(p,1);
x1 = [-1; -1; 1; 1];
x2 = [-1; 1; -1; 1];
X = [x0 x1 x2];

yd = [-1; 1; 1; -1]; % salida deseada

p_cross = 0.35; % proba. de crossover

sigma = 0.45; %tasas de mutaciones para realizaciones distitas

% iteraciones / Generaciones
gen_max = 500; % si nunca converge, corta por aca 
generaciones = zeros(size(vec_individuos));

cant_iteraciones = 100;
idx = 0;

for Nind = vec_individuos
    % itero para distintas tasas de mutacion
    
    

    % Inicializo los pesos de cada individuo. 
    w1_poblacion_inicial = cell(1, Nind); % capa oculta
    w2_poblacion_inicial = cell(1, Nind); % capa salida
    b1_inicial = rand(Nneu, 1) - 0.5; % biasess
    b2_inicial = rand() - 0.5;

    for i = 1:Nind
        w1_poblacion_inicial{i} = rand(Nneu, n) - 0.5;
        w2_poblacion_inicial{i} = rand(Nsal, Nneu) - 0.5;
    end
    
    gen = 0;
    idx = idx + 1;
    display("Cantidad de Individuos: "+ Nind)
    for realizacion = 1:cant_iteraciones
        % entreno la red!
        
        % inicializo los pesos, siempre con el mismo valor random generado
        % al ppio.
        w1_poblacion = w1_poblacion_inicial;
        w2_poblacion = w2_poblacion_inicial;
        % bias inicializo
        b1 = b1_inicial;
        b2 = b2_inicial;
        
        error = 1;
        tolerancia = 0;
        fitness_d = 1 - tolerancia / 4; % fitness target
        fitness_elite = 0;
        fitness_max = fitness_elite;
        
        % calculo inicialmente el fitness:
        for individuo = 1: Nind
            error = 0;

            for mu = 1:p
                x = X(mu, :)';

                h1 = w1_poblacion{individuo} * x + b1;
                V1 = signo(h1);

                h2 = w2_poblacion{individuo} * V1 + b2;
                V2 = signo(h2);

                error = error + (yd(mu) - V2)^2;
            end

            error = error / p;
            vec_fitness(individuo) = 1 - error / 4;
        end

        [fitness_elite, elite_idx] = max(vec_fitness);

        if fitness_elite > vec_fitness(end)
            fitness_max = [fitness_max fitness_elite];
        else
            fitness_elite = fitness_max(end);
            fitness_max = [fitness_max fitness_elite];
        end
        
        % aprendizaje de la red:
        gen = 1;
        while (fitness_elite < fitness_d) && (gen < gen_max)

            % me guardo los indiviudos de elite
            w1_elite = w1_poblacion{elite_idx};
            w2_elite = w2_poblacion{elite_idx};

            p_repr = vec_fitness .* 1 / sum(vec_fitness); % proba de reproduccion

            idx_repr = randsample(1:Nind, Nind, true, p_repr);

            % inicializo la poblacion de reproducidos
            w1_poblacion_repr = cell(1, Nind);
            w2_poblacion_repr = cell(1, Nind);
            for i = 1:Nind
                w1_poblacion_repr{i} = zeros(Nneu, n) - 0.5;
                w2_poblacion_repr{i} = zeros(Nsal, Nneu) - 0.5;
            end

            for individuo1 = 1:(Nind - 1)
                for individuo2 = idx_repr
                    w1_poblacion_repr{individuo1} = w1_poblacion{individuo2};
                    w2_poblacion_repr{individuo1} = w2_poblacion{individuo2};

                    % crossover
                    r = rand();
                    if r < p_cross
                        % se acepta el cruce entre individuos
                        % con uno distinto al individuo actual
                        individuo3 = idx_repr(randi(Nind));
                        while(individuo3 == individuo2)
                            individuo3 = idx_repr(randi(Nind));
                        end
                        % Se hace la cruza:
                        w1_poblacion_repr{individuo1} = w1_poblacion{individuo3};
                        w2_poblacion_repr{individuo1} = w2_poblacion{individuo3};
                    end

                end
            end

            % mutacion
            w1_delta = cell(1, Nind);
            w2_delta = cell(1, Nind);
            for i=1:Nind
                w1_delta{i} = generar_w_mutacion(Nneu, n, sigma);
                w2_delta{i} = generar_w_mutacion(Nsal, Nneu, sigma);
            end    

            for i=1:Nind
                w1_poblacion_repr{i} = w1_poblacion_repr{i} + w1_delta{i};
                w2_poblacion_repr{i} = w2_poblacion_repr{i} + w2_delta{i};
            end
            w1_poblacion_repr{Nind} = w1_elite;
            w2_poblacion_repr{Nind} = w2_elite;

            w1_poblacion = w1_poblacion_repr;
            w2_poblacion = w2_poblacion_repr;


            % recalcula el fitness 
            vec_fitness = zeros(Nind , 1);
            error = 0;

            for individuo = 1: Nind
                error = 0;

                for mu = 1:p
                    x = X(mu, :)';

                    h1 = w1_poblacion{individuo} * x + b1;
                    V1 = signo(h1);

                    h2 = w2_poblacion{individuo} * V1 + b2;
                    V2 = signo(h2);

                    error = error + (yd(mu) - V2)^2;
                end

                error = error / p;
                vec_fitness(individuo) = 1 - error / 4;
            end

            [fitness_elite, elite_idx] = max(vec_fitness);

            if fitness_elite > vec_fitness(end)
                fitness_max = [fitness_max fitness_elite];
            else
                fitness_elite = fitness_max(end);
                fitness_max = [fitness_max fitness_elite];
            end

            gen = gen + 1;
            %display("Generacion: " + gen);
        end
        generaciones(idx) = generaciones(idx) + gen;
    end
    
    generaciones(idx) = generaciones(idx) /  cant_iteraciones;
    % siguiente sigmas
end

plot_d = figure();
plot(vec_individuos, generaciones, 'linewidth', 3, 'color', 'black');
xlabel('Cantidad de Individuos'); ylabel('Generaciones'); 
legend("sigma = " + sigma + newline  + "Proba de crossover = " + p_cross , 'location', 'northwest');
set(gca,'fontsize', 12); 
set(plot_d, 'PaperSize', [20 10]); 
print(plot_d,'Resultados/Ejercicio5/genVSind.pdf','-dpdf');

