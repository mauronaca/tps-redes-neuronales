clc; close all;
imagesPath = "./Imagenes";

%% 45x60
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
    figure()
    imshow( reconstruirImg(hSync, 45, 60));
    figure()
    imshow(reconstruirImg(hAsync, 45, 60));
    display("Error con actualización sincronica:" + sum(hSync - P(:,i))/N);
    display("Error con actualización asincronica:" + sum(hAsync - P(:,i))/N);
end

% Ahora pruebo con la imagen de la paloma con ruido (pixeles negros
% insertados al azar): 
palomaRuido1 = imread('./Imagenes/palomaRuido1.bmp');
P_palomaRuido1 = pixel2Estado(palomaRuido1);
P_paloma = pixel2Estado(paloma);
h = ejecutarRedSync(W, P_palomaRuido1);
figure()
imshow(reconstruirImg(h, 45, 60)); title('Salida con entrada de paloma con ruido');
err = mean(h-P_paloma);
%disp('Error para la paloma con ruido: ' + err);

%% 50x50
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
    figure()
    imshow( reconstruirImg(hSync, 50, 50));
    figure()
    imshow(reconstruirImg(hAsync, 50, 50));
    display("Error con actualización sincronica:" + sum(hSync - P(:,i))/N);
    display("Error con actualización asincronica:" + sum(hAsync - P(:,i))/N);
end

% Ahora la entrada es la imagen del panda con pixeles negros insertados:
pandaRuido1 = imread('./Imagenes/pandaRuido1.bmp');
P_pandaRuido1 = pixel2Estado(pandaRuido1);
P_panda = pixel2Estado(panda);
h = ejecutarRedSync(W, P_pandaRuido1);
figure()
imshow(reconstruirImg(h, 50, 50)); title('Salida con entrada de panda con ruido');
err = mean(h-P_panda);
