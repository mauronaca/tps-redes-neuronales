% Paso los pixeles 0 o 1 a -1 y +1 y ademas hago un reshape
function P = pixel2Estado(img)
    dimension = size(img);
    N = dimension(1) * dimension(2); 
    P = reshape(img, [], 1);
    P = ones(N, 1) - 2 .* P;
end