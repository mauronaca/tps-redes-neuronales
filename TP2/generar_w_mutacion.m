
function w = generar_w_mutacion( Nsig, Nprev, sigma)
    C = sigma * (ones(Nprev, Nprev) - eye(Nprev, Nprev))+eye(Nprev, Nprev);
    w = signo(mvnrnd(zeros(1, Nprev), C, Nsig));
end

% w1 Nneu x n
% w2 Nsal x Nneu