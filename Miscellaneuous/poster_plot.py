
from random import seed
from random import random
seed(1)

interval = np.zeros((20,1))
s = 0
e = 1
for i in range(len(interval)):
    dif = e-s
    if i == 0:
        interval[i] = dif/20
    else:
        interval[i] = interval[i-1] + dif/20

ninterval = np.zeros((20,1))
s = 1
e = 0
for i in range(len(interval)):
    dif = e-s
    if i == 0:
        ninterval[i] = dif/20
    else:
        ninterval[i] = interval[i-1] + dif/20


a1 = 290
ai1 = 20
a2 = 150
ai2 = 20
a3 = 260
ai3 = 20
a4 = 45
ai4 = 20
a5 = 380
ai5 = 20
a6 = 380
ai6 = 20
x1 = np.zeros((a1,1))  # Q
xi1 = np.zeros((20,1)) #+ interval
x2 = np.zeros((a2,1))  #+ 1 + np.random.normal(.1, .2, (a2,1)) # E
xi2 = np.zeros((20,1)) #+ ninterval
x3 = np.zeros((a3,1))  # Q
xi3 = np.zeros((20,1)) + interval
x4 = np.zeros((a4,1))  + 1 + np.random.normal(.1, .2, (a4,1)) # p
xi4 = np.zeros((20,1)) + ninterval#+ 1 + np.random.normal(.1, .2, (20,1)) 
x5 = np.zeros((a5,1))  #+ 1 + np.random.normal(.1, .2, (a5,1)) # t
xi5 = np.zeros((20,1)) #+ ninterval
x6 = np.zeros((a6,1))  # Q
xi6 = np.zeros((20,1))

xs = np.concatenate((x1, xi1, x2, xi2, x3, xi3, x4, xi4, x5, xi5, x6, xi6)) 

book = a1 + a2 + a3 + a4 + a5 + a6 + 120
x = np.random.uniform(0,.5, size=(book,1)) + .5

#x = np.concatenate((x1, x2))

plotx = x + xs
plt.plot(plotx)
plt.grid(False)
figFile = '/Users/marjanfarahbod/Documents/talks/RSGDREAM2022/fig3.pdf'
plt.savefig(figFile)
plt.show()


a1 = 300
a2 = 170
a3 = 280
a4 = 65
a5 = 400
a6 = 400
x1 = np.random.uniform(0,1, size=(a1,1)) 
x2 = np.random.uniform(0,1, size=(a1,1)) 
x3 = np.random.uniform(0,1, size=(a1,1)) 
x4 = np.random.uniform(0,1, size=(a1,1)) 
x5 = np.random.uniform(0,1, size=(a1,1)) 
x6 = np.random.uniform(0,1, size=(a1,1))

