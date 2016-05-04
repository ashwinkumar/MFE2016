import LeastSquareMC as lsMC
import numpy as np
x = np.array([1,2,3,4,5,6])
y = list(lsMC.Hermite(4,x))
print(y)