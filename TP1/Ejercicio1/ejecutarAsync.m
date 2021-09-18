
% Actualizar la red de manera ascincronica
function h = ejecutarAsync(W, entrada)
    dimensionW = size(W);
    N = dimensionW(1);
    h = entrada; % El estado de las neuronas lo pongo inicialmente igual al patron de entrada
    % Para cada neurona -i- la actualizo segun el estado del resto
    for i = 1:N
        % Sumatoria:
        for j = 1:N
            h(i) = h(i) + W(i,j) * sign(h(j));
        end
        h(i) = sign(h(i));
    end
end