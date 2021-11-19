clc; close all; clear all;
imgsRuido = dir('Imagenes/Ruido/*.bmp');

%% Red para imagenes de 45x60
% --------------------------
% Red para set de iamgenes de 45x60
% --------------------------
paloma = imread('./Imagenes/paloma.bmp');
quijote = imread('./Imagenes/quijote.bmp');
torero = imread('./Imagenes/torero.bmp');

set1(:,:,1) = paloma;
set1(:,:,2) = quijote;
set1(:,:,3) = torero;  

[P, W, N, Np] = entrenarRed(set1);

% Compruebo si la red devuelve el mismo patron al inicializarla con ese patron:
disp('Set de imagenes de 45x60:');
for i = 1:Np
    hSync = ejecutarRedSync(W, P(:,i));
    hAsync = ejecutarAsync(W, P(:, i));
    pathName = strcat('Resultados/Ejercicio1/set45x60', num2str(i), '.bmp');
    imwrite(imfuse(reconstruirImg(P(:,i), 45, 60), reconstruirImg(hAsync, 45, 60), 'montage'), pathName);
    display("Error con actualización sincronica:" + sum(hSync - P(:,i))/N);
    display("Error con actualización asincronica:" + sum(hAsync - P(:,i))/N);
end

% Ahora pruebo con ruido (pixeles negros
% insertados al azar): 
for i = 1:6
    % Deberia haber separado cada set en carpetas
    % distintas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if strcmp(imgsRuido(i).name,'perroRuido1.bmp')
        continue
    elseif strcmp(imgsRuido(i).name,'pandaRuido1.bmp')
        continue
    elseif strcmp(imgsRuido(i).name,'vRuido1.bmp')
        continue;
    elseif strcmp(imgsRuido(i).name,'perroRuido.bmp')
        continue;
    end
    display(strcat('Imagenes/Ruido/',imgsRuido(i).name))
    currentImage = imread(strcat('Imagenes/Ruido/',imgsRuido(i).name));
    P = pixel2Estado(double(currentImage));
    hAsync = ejecutarAsync(W, P);
    hSync = ejecutarRedSync(W, P);
    pathName = strcat('Resultados/Ejercicio1/set45x60ruido', num2str(i), '.bmp');
    imwrite(imfuse(currentImage, reconstruirImg(hSync, 45, 60), 'montage'), pathName);
    err = mean(hAsync-P);
    disp("Error con ruido: " + err);
end

% Con la mitad de la imagen del quijote y su otra mitad negro.
quijoteRuido = imread('./Imagenes/Ruido/quijoteRuido.bmp');
P_quijoteRuido = pixel2Estado(quijoteRuido);
h = ejecutarAsync(W, P_quijoteRuido);
pathName = 'Resultados/Ejercicio1/set45x60MitadNegro.bmp';
imwrite(imfuse(quijoteRuido, reconstruirImg(h, 45, 60), 'montage'), pathName)
err = mean(h-pixel2Estado(quijote));
disp("Error con Quijote tapado: " + err);

% Imagen negra:
negro = imread('./Imagenes/Ruido/negro.bmp');
P_negro = pixel2Estado(double(negro));
h = ejecutarAsync(W, P_negro);
pathName = 'Resultados/Ejercicio1/set45x60negro.bmp';
imwrite(imfuse(negro, reconstruirImg(h, 45, 60), 'montage'), pathName)

%% Red para imagenes de 50x50
% --------------------------
% Red para set de imagenes de 50x50
% --------------------------
panda = imread('./Imagenes/panda.bmp');
perro = imread('./Imagenes/perro.bmp');
v = imread('./Imagenes/v.bmp');

set2(:,:,1) = panda;
set2(:,:,2) = perro;
set2(:,:,3) = v;  

[P, W, N, Np] = entrenarRed(set2);

% Compruebo si la red devuelve el mismo patron cuando le in:
disp('Set de imagenes de 50x50:');
for i = 1:Np
    hSync = ejecutarRedSync(W, P(:,i));
    hAsync = ejecutarAsync(W, P(:, i));
    %figure()
    %imshow( reconstruirImg(hSync, 50, 50));
    pathName = strcat('Resultados/Ejercicio1/set50x50', num2str(i), '.bmp');
    imwrite(imfuse(reconstruirImg(P(:,i), 50, 50), reconstruirImg(hAsync, 50, 50), 'montage'), pathName);
    display("Error con actualización sincronica:" + sum(hSync - P(:,i))/N);
    display("Error con actualización asincronica:" + sum(hAsync - P(:,i))/N);
end

% Con la mitad de la imagen del perro y su otra mitad negro.
perroRuido = imread('./Imagenes/Ruido/perroRuido.bmp');
P_perroRuido = pixel2Estado(perroRuido);
P_mixto = sign( 3*(-P(:, 1) + P(:, 2) + P(:, 3)) + 1);
h = ejecutarAsync(W, P_perroRuido);
pathName = strcat('Resultados/Ejercicio1/set50x50mixto', num2str(1), '.bmp');
imwrite(imfuse(imfuse(perroRuido, reconstruirImg(P_mixto , 50, 50), 'montage'), reconstruirImg(h , 50, 50), 'montage'), pathName);
% En este caso hay un estado espurio entre el perro, v y el panda! A prueba
% y error si hago imshow(perro+v-panda) encuentro que es el mismo patron.
% Lo verifico, el estado espurio deberia ser Pmix = -P1 + P2 + P3
disp("Error para el perro tapado y salida (mixto): " + mean(h-pixel2Estado(perro)));

% Pruebo con el v.bmp
vRuido = imread('./Imagenes/Ruido/vRuido.bmp');
P_vRuido = pixel2Estado(vRuido);
h = ejecutarAsync(W, P_vRuido);
pathName = strcat('Resultados/Ejercicio1/set50x50mitadNegro', num2str(1), '.bmp');
imwrite(imfuse(vRuido, reconstruirImg(h , 50, 50), 'montage'), pathName);
disp("Error para el v tapado y salida: " + mean(h-pixel2Estado(v)));

% Ahora la entrada es la imagen con pixeles negros insertados:
for i = 1:10
    if strcmp(imgsRuido(i).name,'quijoteRuido1.bmp')
        continue;
    elseif strcmp(imgsRuido(i).name,'palomaRuido1.bmp')
        continue;
    elseif strcmp(imgsRuido(i).name,'quijoteRuido.bmp')
        continue;
    elseif strcmp(imgsRuido(i).name,'toreroRuido1.bmp')
        continue;
    elseif strcmp(imgsRuido(i).name,'negro.bmp')
        continue;
    end
    display(imgsRuido(i).name)
    currentImage = imread(strcat('Imagenes/Ruido/',imgsRuido(i).name));
    P = pixel2Estado(double(currentImage));
    hAsync = ejecutarAsync(W, P);
    hSync = ejecutarRedSync(W, P);
    pathName = strcat('Resultados/Ejercicio1/set50x50ruido', num2str(i), '.bmp');
    imwrite(imfuse(currentImage, reconstruirImg(hSync, 50, 50), 'montage'), pathName);
    err = mean(hAsync-P);
    disp("Error con ruido: " + err);
end

% otro estado espurio con el panda
pandaSuperRuido = imread('./Imagenes/Ruido/pandaSuperRuido.bmp');
pathName = strcat('Resultados/Ejercicio1/set50x50pandaSuperRuido', num2str(i), '.bmp');
h = ejecutarAsync(W, pixel2Estado(pandaSuperRuido));
imwrite(imfuse(pandaSuperRuido, reconstruirImg(h , 50, 50), 'montage'), pathName);

%% Red para las 6 imagenes

% Rellenar las imagenes para que tengan 3600 pixeles. Segun la
% teoria, al tener 3600  neuronas, en el mejor de los casos la capacidad de
% la red seria de 36 imagenes.
%

%--- Relleno con el primer metodo: 
% Escalos las imagenes a 60x60
nuevaMedida = [60,60];
paloma = imresize(paloma, nuevaMedida);
quijote = imresize(quijote,  nuevaMedida);
torero = imresize(torero,  nuevaMedida);
panda = imresize(panda,  nuevaMedida);
perro = imresize(perro,  nuevaMedida);
v = imresize(v,  nuevaMedida);

set(:,:,1) = paloma; set(:,:,2) = quijote; set(:,:,3) = torero; 
set(:,:,4) = panda; set(:,:,5) = perro; set(:,:,6) = v;     

% Entreno la red:
[P, W, N, Np] = entrenarRed(set);

% Pruebo la red con los mismos patrones de entrada:
for i = 1:Np
    hAsync = ejecutarAsync(W, P(:, i));
    hSync = ejecutarRedSync(W,P(:, i));
    %figure()
    pathName = strcat('Resultados/Ejercicio1/set60x60', num2str(i), '.bmp');
    imwrite(imfuse(set(:,:,i), reconstruirImg(hSync, nuevaMedida(1), nuevaMedida(2)), 'montage'), pathName);
end

% Pruebo con imagenes con ruido
for i = 1:Np
    currentImage = imresize(imread(strcat('Imagenes/Ruido/',imgsRuido(i).name)), nuevaMedida);
    P = pixel2Estado(double(currentImage));
    hAsync = ejecutarAsync(W, P);
    hSync = ejecutarRedSync(W, P);
    pathName = strcat('Resultados/Ejercicio1/set60x60ruido', num2str(i), '.bmp');
    imwrite(imfuse(currentImage, reconstruirImg(hSync, nuevaMedida(1), nuevaMedida(2)), 'montage'), pathName);
end