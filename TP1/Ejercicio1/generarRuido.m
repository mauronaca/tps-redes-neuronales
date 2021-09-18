clc; close all;
paloma = imread('./Imagenes/paloma.bmp');
quijote = imread('./Imagenes/quijote.bmp');
torero = imread('./Imagenes/torero.bmp');
panda = imread('./Imagenes/panda.bmp');
perro = imread('./Imagenes/perro.bmp');
v = imread('./Imagenes/v.bmp');

rng(0);

%% Genero ruido en las imagenes:
% A cada imagen le agrego pixeles negros con probabilidad de 0.2
palomaRuido = insertarRuido(paloma);
pandaRuido = insertarRuido(panda);
figure()
imshow(palomaRuido);
imwrite(palomaRuido, './Imagenes/palomaRuido1.bmp');
figure()
imshow(pandaRuido);
imwrite(pandaRuido, './Imagenes/pandaRuido1.bmp');

%% Tapar mitad de una imagen:
quijoteRuido = llenarMitad(quijote);
figure()
imshow(quijoteRuido);
imwrite(quijoteRuido, './Imagenes/quijoteRuido.bmp');
perroRuido = llenarMitad(perro);
figure()
imshow(perroRuido);
imwrite(perroRuido, './Imagenes/perroRuido');

%%
pxRuido = ones(10000, 1);
for i = 1 : 10000
    if(rand > 0.5)
        pxRuido(i, 1) = 0;
    end
end
%figure()
%hist(pxRuido, [0, 1]);

% Inserta pixeles negro con proba de 0.2
function imgSalida = insertarRuido(imgEntrada)
    dimension = size(imgEntrada);
    alto = dimension(1); 
    ancho = dimension(2);
    imgSalida = imgEntrada;
    
    for i = 1 : alto
        for j = 1 : ancho
            randNum = rand;
            if(randNum > 0.7)
                imgSalida(i, j) = 0; % Pone un px negro en (i,j)
            end
        end
    end
end

% LLena la mitad de la imagen con pixeles negros:
function imgSalida = llenarMitad(imgEntrada)
    dimension = size(imgEntrada);
    alto = dimension(1);
    ancho = dimension(2);
    imgSalida = imgEntrada;
    
    % Creo que se puede hacer matricial
    for i = round(alto/2) : alto
        for j = 1 : ancho
            imgSalida(i, j) = 0;
        end
    end
end
