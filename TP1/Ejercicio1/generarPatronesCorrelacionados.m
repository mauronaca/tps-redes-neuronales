% La funcion recibe cantidad de filas = Neuronas, cantidad de columnas =
% patrones , y un valor de correlacion. Genera una matriz de N x p normal
% multivariada con media 0 y covarianza de acuerdo al valor del parametro
% corr. En realidad var = ro * sigma_i * sigma_j y ro es la correlacion,
% asi que en verdad el parametro que le paso es la covarianza y no la
% correlacion. Devuelve la matriz con 1 o -1.
function P = generarPatronesCorrelacionados(N, p, corr)
    % Genero la matriz de covarianza:
    cov = zeros(p, p);
    for i = 1:p
        for j= 1:p
            if (i == j)
                cov(i, j) = 1;
            else
                cov(i, j) = corr;
            end
        end
    end
    % Genero columnas de variables normales multivariadas.
    P = sign(mvnrnd(zeros(1, p), cov, N));
end