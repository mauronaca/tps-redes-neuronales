function y = func_vecindad(neurona_i, neurona_ganadora, sigma) 
    % recibe pares de coordenadas
    y = exp(  -norm(neurona_i - neurona_ganadora)^2 / (2 * sigma ^ 2) ) ;
end

