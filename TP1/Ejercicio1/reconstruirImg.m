function img = reconstruirImg(P, ancho, alto)
    img = (ones(ancho, alto) - reshape(P, [ancho, alto])) ./2;
end
    