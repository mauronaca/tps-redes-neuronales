% Esto de manera sincrinoica o sea multiplico el patron por la red y
% obtengo los estados todos juntos. La manera ascincronica es ir
% actualizando de a 1 una neurona.
function h = ejecutarRedSync(W, entrada)
    h = signo(W * entrada); %% Cambiar para q sign(0) de 1 o -1
end