clc; close all; 
N = 10000;
p = 3;
PuniformDiscrete = zeros(N, p);
Pnormal = zeros(N, p);
for paso = 1:p
    PuniformDiscrete(:, paso) = 2 .* randi([0, 1], N, 1) - 1; 
    % Armo la matriz de covarianza
    
end
cov = zeros(paso, paso);
    for i = 1:paso
        for j = 1:paso
            if j == i
                cov(i, j) = 1;
            else
                cov(i, j) = 0.9;
            end
        end
    end
    
    Pnormal =  mvnrnd(zeros(1, paso), cov, N); Pnormal = sign(Pnormal);
display(corr(PuniformDiscrete));
display(corr(Pnormal));