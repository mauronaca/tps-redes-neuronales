close all; clear all; clc;
rng('shuffle');

%% Calculo estadistico de la capacidad del percepetron simple

N = 21; % Cantidad de entradas: 20 + Entrada de Bias
patrones = 1:80;
Nap = zeros(size(patrones));
Nap_var = zeros(size(patrones));
iteraciones = 1; % Cantidad de realizaciones q repito el experimiento para obtener una media
Nrep = 100; 
constAprendizaje = 0.99;
iterMaxAprendizaje = 5000;
iter = 1;

for p = patrones
    Nap_means = zeros(1, iteraciones);

    for muestra = 1:iteraciones
        for i = 1:Nrep
            error = 1;
            cantIterAprendizaje = 0; % Cantidad de iteraciones que tarda en aprender.
            y = zeros(p, 1); % Salida del perceptron ante cada patron.
            w = zeros(N, 1); % Pesos

            % Genero patrones al azari con cada x^mu_i = -1 a 1
            % distribucion. La primer entrada la pongo en 1! Que sería el
            % bias.
            % uniforme
            % tamaño de xp: (p, N) 
            xp = 2 * rand(p, N) - 1;
            xp(:, 1) = 1;

            % Genero las salidas para cada patron. Son p patrones 
            yd = generarMuestrasUnifDiscretas(p, 1);

            % Entreno al perceptron:
            while(error > 0) && (cantIterAprendizaje < iterMaxAprendizaje )
                % Elijo un patron al azar
                mu = randi(p); 
                x = transpose(xp(mu, :));

                % Calculo la salida del perceptron 
                y(mu) = signo(transpose(w) * x);

                % Actualizo w:
                w = w + constAprendizaje * (yd(mu) - y(mu)) * x;

                % Calculo el error:
                error = computarError(y, yd);

                cantIterAprendizaje = cantIterAprendizaje + 1;
            end

            %display("El perceptron convergió en " + cantIterAprendizaje + " iteraciones con error: " + error);

           if(error == 0)
               Nap_means(muestra) = Nap_means(muestra) + 1;
           end

        end
    end
    
    Nap(iter) = sum(Nap_means) / iteraciones;
    Nap_var(iter) = var(Nap_means);
    display("Nap para " + p + " patrones es: " + Nap(iter));
    iter = iter + 1;
end

grafico = figure();
%errorbar(patrones, Nap/Nrep, Nap_var)
plot(patrones, Nap / Nrep, 'linewidth', 3, 'color', 'black')
legend("N = " + N);
grid('on');
xlabel("p - # de patrones", 'fontsize', 14);
ylabel("N_{ap} / N_{rep}", 'fontsize', 14);
set(gca,'fontsize', 14);
set(grafico,'PaperSize',[20 10]); %set the paper size to what you want  
print(grafico,'Resultados/Ejercicio2/capacidad','-dpdf') % then print it