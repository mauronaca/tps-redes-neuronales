
function y = batchSizes(max, p)
    iter = 1;
    for i = 1:max
        if rem(p, i) == 0
            y(iter) = i;
            iter = iter + 1;
        end
    end
end