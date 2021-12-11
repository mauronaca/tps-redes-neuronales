function err = computarError(y, yd)
    err = sum((yd - y) .^ 2)/size(yd,1);
end