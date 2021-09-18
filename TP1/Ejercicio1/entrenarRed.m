function [P, W, N, Np] = entrenarRed(set)
    dimension = size(set);
    N = dimension(1) * dimension(2);
    Np = dimension(3);
    P = zeros(N, Np);
    
    for i = 1:Np
        P(:,i) = reshape(set(:,:,i), [], 1);
    end
    
    P = ones(N, Np) - 2 .* P; % Transformo los 0 en +1 y los 1 en -1
    W = P * P' - Np * eye(N);
end
