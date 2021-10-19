# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 17:59:15 2021

@author: mnaca
"""

'''
------------------------------------------------------------------------------
1) Implemente un perceptrón simple que aprenda la función lógica AND de 2 y de 4
entradas. Lo mismo para la función lógica OR. Para el caso de 2 dimensiones, grafique
la recta discriminadora y todos los vectores de entrada de la red.
------------------------------------------------------------------------------
'''

'''
------------------------------------------------
Perceptron simple para funcion AND de 2 entradas
------------------------------------------------
    x1 ｜ x2  ｜ yd
   ----------------
   -1  ｜ -1  ｜ -1
   -1  ｜  1  ｜ -1
    1  ｜ -1  ｜ -1
    1  ｜  1  ｜  1
'''

import numpy as np
from time import time
import matplotlib.pyplot as plt

np.random.seed(seed = int(time()))

def computarError(y, yd):
    p = np.size(yd)
    sum = 0
    for i in range(0, p):
        sum += (yd[:, i] - y[:, i]) ** 2
    return float(sum)

def signo(h):
    return np.sign(2 * h + 1)

p = 4 # Cantidad de patrones
n = 2 # cantidad de entradas

# Patrones. tamaño (4, 2) 4 patrones, 2 entradas
xp = np.array([[-1, -1], [ -1, 1], [1, -1], [1, 1]]).T
# Salida deseada para cada patron
yd = np.array([[-1, -1, -1, 1]]) # Clasificacion, -1 o 1

# Inicio w al azar. Como son 2 entradas w es de tamaño (1,2). w =  w1, w2
#w = np.random.rand(n).reshape(n, 1)
w = np.ones(n).reshape(n ,1)

w0 = np.array([[1]])
error = 1
y = np.zeros(p).reshape(1, p)
constanteDeAprendizaje = 0.99
cantIter = 0

while (error > 0):
    # Elijo un patron al azar
    mu = np.random.randint(0, p) # Elijo un patron al azar
    print('Patron: {}'.format(mu))
    
    # Me quedo con dicho patron y se lo asigno a x
    x = xp[:, 0] 
    x = x.reshape(n, 1) # Reshape queda un vector columna de n filas
    print('Patron actual: {}'.format(x))
    print('w actual: {}'.format(w))
    
    # Calculo la salida del perceptron con la entrada de x. 
    y[0, mu] = signo(np.dot(w.T, x)) # w lo transpongo para que de vector fila
    
    # Actualizo w : aprendizaje con la salida deseda del patron mu y la salida calculada.
    w = w + constanteDeAprendizaje * x * (yd[0, mu] - y[0, mu])
    print('w actualizado: {}'.format(w))
    print('Salida deseada: {}'.format(yd))
    print('Salida del perceptron: {}'.format(y))
    
    # Computo el error entre las salidas deseadas y las salidas del eprceptron para cada patron
    error = computarError(y, yd)
    print('Error: {}'.format(error))
    
    # Contador de iteraciones
    cantIter += 1
    print('\n')

print('Cantidad de iteraciones realizadas: {}'.format(cantIter))

w = np.concatenate((w0, w), axis = 0) # Vector columna

# Plot
x1 = np.linspace(-1, 1, 100).reshape(1, 100)
x2 = - w[0, 0] / w[2, 0] - x1 * w[1, 0] / w[2, 0]


plt.figure(figsize = (10,10))
plt.plot(x1, x2, '-o')
plt.scatter(xp[0,:], xp[1,:], linewidths = 10, color = 'black')
plt.grid(True)
plt.show()