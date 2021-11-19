function y = generarMuestrasUnifDiscretas(N, p)
    rng('shuffle');
    y = zeros(N, p);
    for i = 1:N
        for j = 1:p
            r = rand();
            if(r < 0.5)
               y(i, j) = 1; 
            else
                y(i, j) = -1;
            end  
        end
    end
   
end

