import numpy as np
import mpT3libs

shp3 = [3,3,3]

(El1,El2,El3,Ev1,Ev2,Ev3) = mpT3libs.T3_rnd(shp3)

(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE) = mpT3libs.mpT3_to_T3LE( El1,El2,El3,Ev1,Ev2,Ev3, -1 )

(El1new,El2new,El3new,Ev1new,Ev2new,Ev3new,validMask) = mpT3libs.mpT3LE_to_T3(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE, -1)

print("El1: ")
print(El1)
print("")
print("T(El1): ")
print(T11LE)
print("")
print("invT(T(El1)): ")
print(El1new)
