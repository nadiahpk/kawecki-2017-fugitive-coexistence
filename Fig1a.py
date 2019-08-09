# see if I can replicate Kawecki's (2017) model
# try to verify behaviour of solid line of Figure 1a
# the solid line shows the minimum cost of plasticity that allows the adaptable species to invade
# so show that a point just below the line can't invade and above the line can (and coexistence)

import numpy as np
import matplotlib.pyplot as plt

# patches
# ---

# 101 patches with the optimal trait values ξ_i 
# at regular intervals of 0.02 between -1 and 1

Es = np.arange(-1,1+.02,.02);
no_patches = len(Es)


# parameter values to match from Figure 1a
# ---

s = 0.5     # parameter
v = 0.1     # genetic variaton
m = 0.4     # dispersal probability
p = 0.01    # probability that a patch goes locally extinct (also called "e")

# for the above parameters, c ~ 0.1 in order for adaptable species to invade
c = 0.07    # should not be able to invade here -- yep
c = 0.13    # should be able to invade here -- yep
c = 0.20    # should be able to invade here but will not end in coexistence 
            # -- no, did end in coexistence. Perhaps invasion of plastic not possible though? See Fig1b.py for checking this.


# initialisations
# ---

# To test if the adaptable species can invade, 
# initiated with q_i = 0.001 and x_i = 0 in all patches

qs = np.array( [0.001] * no_patches )
xs = np.array( [0] * no_patches )


# follow evolutionary trajectory
# ---

# The adaptable species was considered to have invaded if its mean relative abundance across
# patches exceeded 0.1 or if the geometric mean of the rate of increase in mean relative abun-
# dance between generations 900 and 1000 was greater than 1. Conversely, competitive exclu-
# sion (i.e., inability to invade) was concluded if became smaller than 10 –6 or if the geometric
# mean of its rate of increase between generations 900 and 1000 was smaller than 1


qb_stor = [ np.mean(qs) ]
for t in range(1000):

    # calculate fitness of each type in each patch
    # ---

    # w_Pi = exp(–c), (3)
    w_P = np.exp(-c)

    # w_Ai ~ f_i (x) = exp[–s(x – ξ)^2 ], (4)
    w_As = np.exp( -s * (xs - Es)**2 )


    # apply catastrophe
    # ---

    # patches that suffer a catastrophe marked true
    catast = np.random.random_sample(no_patches) <= p
    prop_surv = ( no_patches - sum(catast) ) / no_patches


    # mean genotype after selection but before dispersal
    # ---

    xms = np.array([ x + 2*s*v*(E-x) if not c else 0 for x, E, c in zip(xs, Es, catast) ])


    # dispersal and reproduction
    # ---

    # q_i^* = qs*w_As / ( qs*w_As + (1-qs)*w_Ps )
    qms = np.array([ q*w_A / (q*w_A + (1-q)*w_P) if not c else 0 for q, w_A, c in zip(qs, w_As, catast) ])

    # \bar{q}^*
    sum_qms = sum(qms)
    qb = sum_qms / no_patches

    # \bar{x}^*
    xb = sum(qms*xms) / sum_qms

    '''
    # taking his equations literally

    # Equation 7 xd = ( (1-m)*qms*xms + m*(1-p)*qb*xb ) / ( (1-m)*qms + m*(1-p)*qb )
    xds = np.array([ ( (1-m)*qm*xm + m*(1-p)*qb*xb ) / ( (1-m)*qm + m*(1-p)*qb ) if not c else xb for qm, xm, c in zip(qms, xms, catast) ])

    # Equation 2 qd = ( (1-m)*qms + m*(1-p)*qb ) / ( 1 - m + m*(1-p) )
    qds = np.array([ ( (1-m)*qm + m*(1-p)*qb ) / ( 1 - m + m*(1-p) ) if not c else qb for qm, c in zip(qms, catast) ])
    '''

    # my own variant, bc we know the proportion who survived

    # Equation 7 xd = ( (1-m)*qms*xms + m*(1-p)*qb*xb ) / ( (1-m)*qms + m*(1-p)*qb )
    xds = np.array([ ( (1-m)*qm*xm + m*prop_surv*qb*xb ) / ( (1-m)*qm + m*prop_surv*qb ) if not c else xb for qm, xm, c in zip(qms, xms, catast) ])

    # Equation 2 qd = ( (1-m)*qms + m*(1-p)*qb ) / ( 1 - m + m*(1-p) )
    qds = np.array([ ( (1-m)*qm + m*prop_surv*qb ) / ( 1 - m + m*prop_surv ) if not c else qb for qm, c in zip(qms, catast) ])

    # update for next iteration
    xs = xds
    qs = qds

    # storage
    qb_stor.append(qb)


# plot the mean proportion of adaptable species with time
# ---

plt.plot(qb_stor)
plt.xlabel('iterations')
plt.ylabel(r'mean propn adaptable species, $\bar{q}^*$')
plt.title(r'Fig 1a with $e=' + str(p) + '$ at point $m=' + str(m) + ', c=' + str(c) +'$')
plt.tight_layout()
plt.savefig('Fig1a_c' + str(int(np.floor(c*100))) + '.png')
plt.close()
