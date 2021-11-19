function y = computarSalidaBias(x, w, b)
    n = size(w, 2); % #capas
    V{1} = x;
    %N = zeros(1, n); % #neuronas p/capa
    for m = 1:n
        N(m) = size(w{m}, 1);
        h{m} = zeros(N(m), 1);
        V{m+1} = tanh(h{m});
    end
  
    for m = 1:n
        for i = 1:N(m)
            h{m}(i) = w{m}(i, :) * V{m} + b{m}(i);
            V{m+1}(i) = tanh(h{m}(i));  
        end
    end
    
    y = V{n+1};
end