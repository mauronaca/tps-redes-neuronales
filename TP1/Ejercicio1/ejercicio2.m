close all; clc; clear all;

%% Compruebo estadisticamente la capacidad de la red
% Genero p patrones aleatorios de N neuronas, luego entreno una red y
% computo el error entre el patron y la salida de la red ante el mismo
% patron. Si el error es menor al umbral p_error entonces repito los
% primeros pasos aumentando la cantidad de patrones. La cantidad de patrones que resultan cuando el error supere
% el umbral, es la capacidad.

listaErrores = [0.001, 0.0036, 0.01, 0.05, 0.1];
capacidadMedia = zeros(1, 5);
capacidadVar = zeros(1, 5);
cant_iteraciones = 10;
N = 100;

for i = 1:5
    p_error = listaErrores(i);
    p = 1;
    error = 0;
    c = zeros(1, cant_iteraciones);
    
    % Repito el expermiento para sacar media y varianza.
    for muestra = 1:cant_iteraciones
        p = 1;
        error = 0;
        while (error < p_error)
            % Genero los patrones pseudoaleatorios. 
            P = 2 .* randi([0, 1], N, p) - 1;

            % Entreno la red, calculo la matriz de pesos:
            W = P * P' - p * eye(N);

            % Ejecuto la red y calculo el error:
            error = mean(mean((sign(W * P) - P) ~=0 ));
            p = p + 1;
        end
        c(1, muestra) = (p - 1)/N;
        %disp("Capacidad para la iteracion " + muestra + " = " + c(1, muestra));
    end
    
    capacidadMedia(i) = mean(c);
    capacidadVar(i) = sqrt(var(c));
    disp("Capacidad total con P_error " + p_error + " = " + capacidadMedia(i) + " con varianza " + capacidadVar(i));
    
end
h = figure()
errorbar(listaErrores, capacidadMedia, capacidadVar, 'linewidth', 2, 'color', 'red');
title("Capacidades en función de la probabilidad de error para "+N+" neuronas");
legend('Capacidad');
grid('on');
xlabel('Probabilidad de error admitida');
ylabel('Capacidad [p_{max} / N]');
set(gca,'fontsize', 14);
saveas(gcf, 'Resultados/Ejercicio2/capacidad.png');
set(h,'PaperSize',[20 10]); %set the paper size to what you want  
print(h,'Resultados/Ejercicio2/grafico1','-dpdf') % then print it

%hist(P, [-1,1])

%% Con distintas correlaciones
% corr(P) calcula la correlacion para cada par de columnas
% El metodo para agregarle correlacion a los patrones es: generar matrices
% normales multivariadas con una matriz de covarianza que depende de la
% correlacion que quiera darle. La diagonal de la matriz de cov será 1 y
% los demas elementos iran desde 0 hasta 0.8. Luego de generar los patrones
% calculo la correlacion. 
% si paso = p  => correlacion es 100% 

N = 100;
cant_pasos = 50;
img_idx = 1;
cantIteraciones = 200;
for probaError = listaErrores
    
    capacidad = zeros(1, cantIteraciones);
    capacidades = zeros(1, cant_pasos);
    idx = 1;
    correlaciones = linspace(0.01, 0.99, cant_pasos);
    %correlaciones = [linspace(0,0.5,50), linspace(0.51,0.8,50)];
    correlacion = zeros(1, cantIteraciones);
    varianza = zeros(1, cant_pasos);
    for ro = correlaciones
        % Este expermiento habria que repetirlo varias veces y sacar media y
        % varianza.
        for i = 1:cantIteraciones
            p = 1;
            error = 0;
            while(error < probaError)
                % Genero patrones correlacionados:
                P = generarPatronesCorrelacionados(N, p, ro);

                % Entreno la red, calculo la matriz de pesos:
                W = P * P' - p * eye(N);

                %Actualizo con la verdadera correlacion q da
                correlacion(i) = mean(mean(corr(P)));

                % Ejecuto la red y calculo el error:
                error = mean(mean((sign(W * P) - P) ~=0 ));
                p = p + 1;
            end  
            capacidad(i) = (p - 1) / N;
        end
        correlaciones(idx) = mean(correlacion);
        capacidades(idx) = mean(capacidad);
        varianza(idx) = sqrt(var(capacidad));
        idx = idx + 1;
    end
    h = figure()
    errorbar(correlaciones, capacidades, varianza , 'linewidth', 2, 'color', 'red');
    title("Capacidades en función de la correlación para "+N+" neuronas y proba. de error admitida "+probaError);
    legend('Capacidad');
    grid('on');
    xlabel('Correlación entre los patrones enseñados');
    ylabel('Capacidad [p_{max} / N]');
    set(gca,'fontsize', 14);
    pathName = strcat('Resultados/Ejercicio2/capacidad', num2str(img_idx));
    %saveas(gcf, 'Resultados/Ejercicio2/capacidad.png');
    set(h,'PaperSize',[20 10]); %set the paper size to what you want  
    print(h,pathName,'-dpdf') % then print it
    img_idx = img_idx + 1;
end

