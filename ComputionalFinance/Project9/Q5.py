import numpy as np
import InterestRateModel as IRM
import  PrePaymentModel as PPM
import timeit
import matplotlib.pyplot as plt
from tabulate import tabulate

T = 30
WAC = 0.08
PV0 = 100000

### CIR Model
r0 = 0.078
kappa = 0.6
rbar_array = np.arange(0.03,0.099,0.01)
sigma = 0.12
dT = 1/360
nsims = 1000

IO = np.zeros(len(rbar_array))
PO = np.zeros(len(rbar_array))
#### Question 5
'''
Consider the MBS described above and the IO and PO tranches.
Use the Numerix-Prepayment Model and price the IO and PO tranches for: 𝑟̅ in 0.03 to 0.09 range, in increments of 0.01.
'''

start_time1 = timeit.default_timer()
for i in range(len(rbar_array)):
    cir = IRM.CIR(r0= r0,r_bar=rbar_array[i],T=T+10, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
    r_sim1 = cir.EulerDiscretize_R()
    numerixPPM = PPM.NumerixPrepayment(PV0, WAC, r_sim1, T)
    numerixPPM.setInterestRateModel(cir)
    IO[i] = numerixPPM.MBS_IOPrice()
    PO[i] = numerixPPM.MBS_POPrice()
    del numerixPPM
    del cir

result = np.zeros([3,len(rbar_array)])
result[0,:] = rbar_array
result[1,:] = IO
result[2,:] = PO

ax = plt.subplot(1, 1, 1)

plt.xlabel("rbar")
plt.ylabel("MBS Tranch Prices")
plt.title("IO/PO Price vs rbar of CIR Model")
plt.plot(rbar_array,IO, label='IO')
plt.plot(rbar_array,PO, label='PO')
ax.legend(loc="center left", bbox_to_anchor=[1, 0.5], ncol=1, shadow=True, prop={'size':6}, title="Leg", fancybox=True)

ax.get_legend().get_title().set_color("red")


plt.savefig("5fig")


elapsed1 = timeit.default_timer() - start_time1
print('Table for IO and PO Prices ****[%f sec]'%elapsed1)

print(tabulate(np.transpose(result), headers=["r","IO","PO"]))
