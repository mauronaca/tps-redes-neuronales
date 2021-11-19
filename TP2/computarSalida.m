function y = computarSalida(x, w, b)
    n = size(w, 2); % #capas
    V{1} = x;
    %N = zeros(1, n); % #neuronas p/capa
    for m = 1:n
        N(m) = size(w{m}, 1);
        h{m} = zeros(N(m), 1);
        V{m+1} = tanh(h{m});
    end
  
    for m = 1:n
        h{m} = w{m} * V{m} + b{m};
        V{m+1} = tanh(h{m});
    end
    
    y = V{n+1};
end