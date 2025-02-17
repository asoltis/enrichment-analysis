import stattests
import time
import numpy as np

pvs = stattests.fishers_exact_test(12,5,29,2,TTmethod='smallP')
print(pvs)
pvs = stattests.fishers_exact_test(255,186,71,183)
print(pvs)

pvs = stattests.hypgeo_test(695,14000,2830,2724)
print(pvs)

# Check timing for 1000 simulations
t0 = time.time()
#for i in range(0,1000):
#    pvs = stattests.fishers_exact_test(1200,1200,1500,1000,TTmethod='smallP')
#print 'Small p method: %f' %(time.clock() - t0)

t0 = time.time()
for i in range(0,1000):
    pvs = stattests.fishers_exact_test(1200,1200,1500,1000,TTmethod='double')
print('Doubling method: %f' %(time.time() - t0))

# FDR
pvals = [1e-3, 1e-5, 1e-10, 0.4, 0.01, 0.25]
pvals = np.array(pvals)
fdrs = stattests.fdr_BH(pvals)
for (p, f) in zip(pvals, fdrs):
    print(p, f)

# Binomial test
prob = stattests.binomial_prob(500, 100, 0.5)
print('Binomail test prob: %g'%(prob))

