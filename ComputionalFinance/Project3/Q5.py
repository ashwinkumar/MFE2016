import myHalton as halton
import myrandgen as rdm
import matplotlib.pyplot as plt
import math
import mywrapper as wp
import statistics as st

def myfunc(x,y):
    done = object()
    n1 = next(x, done)
    n2 = next(y, done)
    while (n1 is not done) and (n2 is not done):
        yield integral(n1,n2)
        n1 = next(x, done)
        n2 = next(y, done)


def integral(x,y):
   d1 = math.exp(-x*y)
   s = math.sin(6*math.pi*x)
   c = wp.mycuberoot(math.cos(2*math.pi*y))
   return d1*(s + c)

n = 100
ugen1 = rdm.LGMRandomGenerator(14)
ugen2 = rdm.LGMRandomGenerator(18)
urand1 =  list(ugen1.generateRdmNumberByGenerator(n))
urand2 =  list(ugen2.generateRdmNumberByGenerator(n))


plt.scatter(urand1, urand2)
plt.ylabel("Y")
plt.xlabel("X")
plt.title("Uniform Random X-Y")
plt.show()


lds1 = halton.Halton(2)
haltonsq1 = list(lds1.generateRdmNumberSequences(n))
lds2= halton.Halton(7)
haltonsq2 = list(lds2.generateRdmNumberSequences(n))


plt.scatter(haltonsq1, haltonsq2)
plt.ylabel("Y")
plt.xlabel("X")
plt.title("Haton(2,7) Sequences X-Y")
plt.show()

lds3 = halton.Halton(4)
haltonsq3 = list(lds3.generateRdmNumberSequences(n))

plt.scatter(haltonsq1, haltonsq3)
plt.ylabel("Y")
plt.xlabel("X")
plt.title("Halton(2,4) Sequences X-Y")
plt.show()

## Integral
n = 10000
lds_2 = halton.Halton(2)
lds_4 = halton.Halton(4)
lds_5 = halton.Halton(5)
lds_7 = halton.Halton(7)


ldseq_2 = lds_2.generateRdmNumberSequences(n)

ldseq_4 = lds_4.generateRdmNumberSequences(n)


ldseq_5 = lds_5.generateRdmNumberSequences(n)

ldseq_7 = lds_7.generateRdmNumberSequences(n)

a1 = list(myfunc(ldseq_2,ldseq_4))
a1_mean = st.mean(a1)
a1_var = st.variance(a1)

ldseq_2 = lds_2.generateRdmNumberSequences(n)

a2 = list(myfunc(ldseq_2,ldseq_7))
a2_mean = st.mean(a2)
a2_var = st.variance(a2)


ldseq_7 = lds_7.generateRdmNumberSequences(n)
a3 = list(myfunc(ldseq_5,ldseq_7))
a3_mean = st.mean(a3)
a3_var = st.mean(a3)

print('Q5.  e) The integral of the function using MC using Halton sequences: ')
print('        a) H(2,4) = %f with var = %f' % (a1_mean,a1_var))
print('        b) H(2,7) = %f with var = %f' % (a2_mean, a2_var))
print('        a) H(5,7) = %f with var = %f' % (a3_mean,a3_var))

ugen_1 = rdm.LGMRandomGenerator(14)
ugen_2 = rdm.LGMRandomGenerator(18)
urand_1 =  ugen_1.generateRdmNumberByGenerator(n)
urand_2 =  ugen_2.generateRdmNumberByGenerator(n)
a4 = list(myfunc(urand_1,urand_2))
a4_mean = st.mean(a4)
a4_var = st.variance(a4)
print('Q5.  e) The integral of the function using MC using Random sequences = %f with var = %f' % (a4_mean,a4_var))


