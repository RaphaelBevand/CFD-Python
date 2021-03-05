import matplotlib.pyplot as plt 
import numpy as np 

x = np.arange(0, 2 * np.pi, 0.1) 
X,Y = np.meshgrid(x,x) 
f1 = np.sin(X) + np.sin(Y) 
f2 = np.cos(X) + np.cos(Y) 

plt.figure() 
C = plt.contourf(f1) 
plt.show() 
for coll in C.collections: 
    plt.gca().collections.remove(coll) 
C = plt.contourf(f2) 
plt.draw() 