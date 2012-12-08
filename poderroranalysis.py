from matplotlib import rc
import os
import matplotlib.pyplot as plt

N = [20, 1] 
dt = 0.01
for i in range(2,4):
    s = "./ns womersley%dd ipcs dt=%g N=%d compare_POD=True" %(i,dt,N[i-2])
    os.system(s)

    errorlist = []
    f = open("POD_error.txt","r")

    for line in f:
        errorlist.append(float(line))

    rc('text', usetex = True)
    rc('font', family='serif')
    plt.plot(range(1,len(errorlist)+1), errorlist, label=r'$E_k$')
    plt.title("Error of " + r'$u_h - \sum_{i=1}^k a_i \phi_i$' + " %dD" %i)
    plt.legend(bbox_to_anchor=(0.9, 0.9), loc=2, borderaxespad=0.)
    plt.xlabel("POD approximation of order k")
    plt.show()
