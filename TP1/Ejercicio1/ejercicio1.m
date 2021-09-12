clc; close all;
imagesPath = "./Imagenes";
%%
% --------------------------
% Red para set de iamgenes de 60x45
% --------------------------

paloma = imread('./Imagenes/paloma.bmp');
quijote = imread('./Imagenes/quijote.bmp');
torero = imread('./Imagenes/torero.bmp');

set1(:,:,1) = paloma;
set1(:,:,2) = quijote;
set1(:,:,3) = torero;  

[P, W, N, Np] = entrenarRed(set1);

% Compruebo si la red devuelve el mismo patron cuando le in:
for i = 1:Np
    h = ejecutarRed(W, P(:,i));
    figure()
    imshow( (ones(45, 60) - reshape(h, [45,60])) ./2 );
    %error = sum(h - P(:,i))/N;
end


%%
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
for i = 1:Np
    h = ejecutarRed(W, P(:,i));
    figure()
    imshow( (ones(50, 50) - reshape(h, [50,50])) ./2 );
    %error = sum(h - P(:,i))/N;
end

%%

function [P, W, N, Np] = entrenarRed(set)
    dimension = size(set);
    N = dimension(1) * dimension(2);
    Np = dimension(3);
    P = zeros(N, Np);
    
    for i = 1:Np
        P(:,i) = reshape(set(:,:,i), [], 1);
    end
    
    P = ones(N, Np) - 2 .* P; % Transformo los 0 en +1 y los 1 en -1
    W = P*P' - Np * eye(N);
end

function h = ejecutarRed(W, entrada)
    h = sign(W * entrada);
end

