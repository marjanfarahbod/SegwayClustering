import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0.01,0.05,100)

y1 = x**(.8)
y2 = x**(.1)
y3 = x**(.01)
y4 = x**(.001)

fig = plt.figure()
plt.plot(x, y1)
plt.plot(x, y2)
plt.plot(x, y3)
plt.plot(x, y4)
plt.show()


fig, axs = plt.subplots(2,2)
axs[0,0].plot(x,y1)
axs[0,1].plot(x,y2)
axs[1,0].plot(x,y3)
axs[1,1].plot(x,y4)

fig.show()
