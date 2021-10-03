clc; clear all; close all;

%% Graficar el error en funcion del porcentaje de conexiones eliminadas

p = 15; % Cantidad de patrones
N = 200; % # de neuronas
cantIteraciones = 100;
conexionesEliminadas = linspace(0, 100, 100); % Porcentajes
cant = size(conexionesEliminadas); cant = cant(2); % Cantidad de porcentajes a borrar
errores = zeros(1, cant);
varianza = zeros(1, cant);
idx = 1;

for porcentaje = conexionesEliminadas
    error = zeros(1, cantIteraciones);
    for iter = 1:cantIteraciones
        % Genero p patrones pseudoaleatorios:
        P = 2 .* randi([0, 1], N, p) - 1;

        % Entreno la red
        W = P * P' - p * eye(N);

        % Elimino una x cantidad de conexiones:
        % Hacerlo al azar preguntando si no estaba en 0.
        x = porcentaje * N * N / 100;
        W = reshape(W, [], 1);
        for i = 1:x
            W(i) = 0;
        end
        W = reshape(W, N, N);
        error(iter) =  mean(mean((signo(W * P) - P) ~=0 ));
    end
    
    errores(1, idx) = mean(error);
    varianza(1, idx) = sqrt(var(error));
    %display("Error = " + errores(1, idx) );
    %display("Varianza = " + sqrt(var(error)) );
    idx = idx + 1;
end

figure()
errorbar(conexionesEliminadas, errores, varianza, 'color', 'red', 'linewidth', 2);
title("Proba. de error en función del porcentaje de conexiones eliminidas para "+ N +" neuronas y "+ p +" patrones");
legend("Error para " + p + " patrones");
grid('on');
xlabel('Porcentaje de conexiones entre neuronas eliminadas');
xtickformat('percentage');
ylabel('Probabilidad de error');
set(gca,'fontsize', 14);
saveas(gcf, 'Resultados/Ejercicio3/erroresPercent.png');


%% Graficar la capacidad de la red en funcion del porcentaje de neuronas eliminadas:

p = 1; % Cantidad de patrones
N = 150;
cantIteraciones = 100;
conexionesEliminadas = linspace(0, 100, 100); % Porcentajes
cant = size(conexionesEliminadas); cant = cant(2); % Cantidad de porcentajes a eliminar
errores = zeros(1, cant);
varianza = zeros(1, cant);
capacidades = zeros(1, cant);
idx = 1;
pError = 0.001;

for porcentaje = conexionesEliminadas
    error = zeros(1, cantIteraciones);
    capacidad = zeros(1, cantIteraciones);
    
    % Repito cantIteraciones veces el experimento de calcular la capacidad
    % eliminando un porcentaje de conexiones
    x = porcentaje * N * N / 100;
    for iter = 1:cantIteraciones
        p = 1;
        
        while (error(iter) < pError)
            % Genero p patrones pseudoaleatorios:
            P = 2 .* randi([0, 1], N, p) - 1;

            % Entreno la red
            W = P * P' - p * eye(N);

            % Elimino una x cantidad de conexiones:
            W = reshape(W, [], 1);
            for i = 1:x
                W(i) = 0;
            end
            W = reshape(W, N, N);
            error(iter) =  mean(mean((signo(W * P) - P) ~=0 ));
            p = p + 1;
        end
        
        capacidad(iter) = (p - 1) / N;
    end
    
    errores(1, idx) = mean(error);
    varianza(1, idx) = sqrt(var(capacidad));
    capacidades(1, idx) = mean(capacidad);
    %display("Error = " + errores(1, idx) );
    %display("Varianza = " + sqrt(var(error)) );
    %display("Capacidad = " + capacidades(1, idx));
    idx = idx + 1;
end

figure()
errorbar(conexionesEliminadas, capacidades, varianza, 'linewidth', 2, 'color', 'red');
title("Capacidad en función del porcentaje de conexiones eliminidas para "+N);
grid('on');
xlabel('Porcentaje de conexiones entre neuronas eliminadas');
xtickformat('percentage');
ylabel('Capacidad [p_{max} / N]');
set(gca,'fontsize', 14);
saveas(gcf, 'Resultados/Ejercicio3/capacidadPercent.png');
