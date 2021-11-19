# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 18:42:10 2021

@author: mnaca
"""

import numpy as np

def redondear(num, digitos):
    # Recibe un numero entre 1 y 0 y una cantidad de digitos. Redondea para arriba o 
    # para abajo.
    if num >= 1 or num <= 0:
        raise ValueError('Numero mayor a 1 o menor a 0')
    
    str_num = list(str(num))[2:2 + digitos + 1]
    corte = int(str_num[-1])
    
    if corte > 5:
        str_num[-2] = str(int(str_num[-2]) + 1)
    
    str_num = str_num[0:len(str_num) - 1]
    
    return int(''.join(str_num)) / 10**len(str_num)

def cortar(numero, digitos):
    str_num = str(numero)
    decimal = str_num.split('.')[1]
    entera = str_num.split('.')[0]
    
    num_redondeado = redondear(numero / 10 ** len(list(entera)), digitos) * 10 ** len(list(entera))
    return num_redondeado

assert 0.5 == redondear(0.54, 1)
assert 0.6 == redondear(0.56, 1)
assert 0.55 == redondear(0.554, 2)
assert 0.56 == redondear(0.556, 2)
assert 0.556 == redondear(0.5556, 3)
assert 0.555 == redondear(0.5554, 3)
assert 12.9 == cortar(12.89213123, 3)

