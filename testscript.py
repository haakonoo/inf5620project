from matplotlib import rc
import os
import matplotlib.pyplot as plt


n = 6
N = [2**i for i in range(1,n+1)]
dt = [2**(-i) for i in range(1,n+1)]

path = os.getcwd()
path = os.path.join(path,"error.log")
if(os.path.exists(path)):
    os.remove(path)

for i in range(n):
    s = "./ns test2d ipcs dt=%g N=%d" %(dt[i],N[i])
    os.system(s)

errorlist = []
f = open("error.log","r")
i = 1
for line in f:
    errorlist.append(float(line.split("=")[-1]))
    errorlist[-1] = errorlist[-1]*2**(i)
    i+= 1

rc('text', usetex = True)
rc('font', family='serif')
plt.plot(errorlist, label=r'$E2^i$')
plt.title("Error measured as " + r'$ E = ||u-u_h||_L{^2}([0,T]\times L^2(\Omega))$')
plt.legend(bbox_to_anchor=(0.9, 0.9), loc=2, borderaxespad=0.)
plt.xlabel("Number of iterations i")
plt.show()
