clc;
imagesPath = "./Imagenes";

%%
% --------------------------
% Red para imagenes de 60x45
% --------------------------
paloma = imread('./Imagenes/paloma.bmp');
quijote = imread('./Imagenes/quijote.bmp');
torero = imread('./Imagenes/torero.bmp');

N = 2700;
Np =  3;

% Matriz de patrones
P(:, 1) = reshape(paloma, [], 1);
P(:, 2) = reshape(quijote, [], 1);
P(:, 3) = reshape(torero, [], 1);
%arrayfun(@(x) 1 - 2 * x, P); % Los ceros van a +1 y los unos a -1
P = ones(N, Np) - 2 .* P;
W = P*P' - Np * eye(N);

% Ejecucion:

% Compruebo si la red devuelve el mismo patron:
for i=1:Np
    h = sign(W*P(:, i));
    figure()
    imshow( (ones(45, 60) - reshape(P(:,i), [45,60])) ./2 );
    error = sum(h - P(:,i))/N;
end
%%

%%
% --------------------------
% Red para imagenes de 50x50
% --------------------------
perro = imread('./Imagenes/panda.bmp');
panda = imread('./Imagenes/perro.bmp');
v = imread('./Imagenes/v.bmp');

N = 2500;
Np =  3;

% Matriz de patrones
P(:, 1) = reshape(perro, [], 1);
P(:, 2) = reshape(panda, [], 1);
P(:, 3) = reshape(v, [], 1);
P = ones(N, Np) - 2 .* P;
W = P*P' - Np * eye(N);

% Ejecucion:

% Compruebo si la red devuelve el mismo patron:
for i=1:Np
    h = sign(W*P(:, i));
    figure()
    imshow( (ones(50, 50) - reshape(P(:,i), [50,50])) ./2 );
    error = sum(h - P(:,i))/N;
end
%%