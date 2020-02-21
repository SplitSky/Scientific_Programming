import numpy as np
import matplotlib.pyplot as plt

m = 2
c = 47
x = [1,2,3,4,5,6,7]
x = np.array(x)
y = m*x + c
fit = np.polyfit(x,y,1)
print(fit)
plt.plot(x,y)
plt.show()