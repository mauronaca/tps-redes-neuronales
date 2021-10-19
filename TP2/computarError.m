
function err = computarError(y, yd)
    err = sum((yd - y) .^ 2);
end