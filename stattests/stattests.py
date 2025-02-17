# Python function for computing various statistical tests
# 
# By Anthony Soltis
# 04/07/2017
#
# UPDATES
# 04/24/2018 - included negative binomial probability function (NB_prob)
# 04/16/2019 - included Benjamini-Hochberg FDR function

# Import statements
from __future__ import division
from math import factorial, log, exp, pi
import numpy as np

##################################
# MULTIPLE HYPOTHESIS CORRECTION #
##################################
def fdr_BH(pvals):
    '''
    Compute Benjamini-Hochberg FDR q-values on array of p-values.
    Function returns FDR values in original sort order of input. 
    '''
    total = pvals.size
    si = np.argsort(pvals) 
    pvals_sort = pvals[si]
    
    fdrs0 = total*pvals_sort/range(1,total+1)
    fdrs = []
    # Preserve monotonicity
    for i in range(0,total):
        fdrs.append(min(fdrs0[i:]))
    
    # preserve original sort order 
    pairs = [(i,j) for i, j in enumerate(si)]
    fdrs_resort = np.zeros(pvals.shape)
    for pair in pairs: fdrs_resort[pair[1]] = fdrs[pair[0]]
    return fdrs_resort

####################
# COUNT STATISTICS #
####################

# Fisher's exact test
def fishers_exact_test(a,b,c,d,TTmethod='double'):
    '''
    Compute right, left, and two-tail Fisher's exact test p-values
    according to 2x2 contingency table values:

    -----------------------------------------------
    |         | Outcome 1 | Outcome 2 |    Sum    |
    -----------------------------------------------
    | Group 1 |     a     |     b     |   a + b   |
    -----------------------------------------------
    | Group 2 |     c     |     d     |   c + d   |
    -----------------------------------------------
    |   Sum   |   a + c   |   b + d   | a+b+c+d=n |
    -----------------------------------------------

    with probability:

         /       \ /       \
         | a + b | | c + d |
         |   a   | |   c   |
         \       / \       /
    p = --------------------
              /       \
              |   n   |
              | a + c |
              \       /
    
    INPUTS:
        contingency table values a, b, c, and d as integers
        TTmethod - method for computing two-tail test. Options are:
                   'double' - Simply double minimum of right and left-tails (faster, default). 
                   'smallP' - all p <= p(contingency table) are summed (slower).
    
    OUTPUT:
         3-element list containing left-tail, right-tail, and two-tail p-values.
    '''

    x = a
    M = a+b+c+d
    K = a+b
    N = a+c
    
    return hypgeo_test(x,M,K,N,TTmethod=TTmethod)

# Hypergeometric test
def hypgeo_test(x,M,K,N,TTmethod='double'):
    ''' 
    Compute right, left, and two-tail hypergeometric p-values.
    For right-tail:

                 /   \/       \
        min(K,N) | K || M - K |
          ---    | i || N - i |
          \      \   /\       /
    pv =  /      --------------
          ---        /   \
          i=x        | M |
                     | N |
                     \   /

    where:
        x - number of positives in foreground set
        M - total size of set, or background size
        K - total positives in entire set
        N - total selected in test

    Inputs:
        x, M, K, N as integers
        TTmethod - method for computing two-tail test. Options are:
                   'double' - Simply double minimum of right and left-tails (faster, default). 
                   'smallP' - all p <= p(contingency table) are summed (slower).
        
        Outputs:
        3-element list containing left-tail, right-tail, and two-tail p-values.
    '''         
    # Compute right and left tails, make sure value does not exceed 1
    hg_rt,hg_lt = 0,0
    for i in range(x,min(K,N)+1):
        hg = hypge_pdf(i,M,K,N)
        hg_rt += hg
        if i > x:
            hg_lt += hg
    if hg_rt > 1: hg_rt = 1
    if hg_lt > 1: hg_lt = 1
    hg_lt = 1-hg_lt
    
    # Get two-tail
    # Double one-tail method
    if TTmethod == 'double':
        if hg_lt < hg_rt: two_tail = 2*hg_lt
        elif hg_rt < hg_lt: two_tail = 2*hg_rt

    # Method of summing small p
    elif TTmethod == 'smallP':
        hgsum = 0
        pcrit = hypge_pdf(x,M,K,N)
        for i in range(0,min(K,N)+1):
            try:
                hg = hypge_pdf(i,M,K,N)
                if hg <= pcrit: hgsum+=hg
            except: continue
        two_tail = hgsum

    # Return values
    return [hg_lt,hg_rt,two_tail]

def NB_prob(k,r,p):
    '''
    Compute negative binomial probability with parameters k, r, p according to:

    
        /           \ 
    p = | k + r - 1 | * (1-p)^r * p^k
        |     r     |
        \           /

    '''

    logP = r*log(1-p) + k*log(p)
    xCy = log_xCy(k+r-1,r)

    prob = exp(logP + xCy)
    return prob

####################
# HELPER FUNCTIONS #
####################

def hypge_pdf(x,M,K,N):
    '''
    Compute hypergeomentric probability for given input quantities.
    '''
    # Compute logarithms of "x choose y" values
    a = log_xCy(K,x)
    b = log_xCy(M-K,N-x)
    c = log_xCy(M,N)

    # Exponentiate to put back on normal scale
    hpdf = exp(a + b -c)
    return hpdf

def log_xCy(x,y):
    '''
    Compute the logarithm of "x choose y".
    '''
    # x-value
    if x <= 100:
        fx = log(factorial(x))
    else:
        fx = approx_logn_factorial(x)

    # y-value
    if y <= 100:
        fy = log(factorial(y))
    else:
        fy = approx_logn_factorial(y)

    # x - y value
    if (x-y) <= 100:
        fxy = log(factorial(x-y))
    else:
        fxy = approx_logn_factorial(x-y)

    # Assure values are non-negative
    fx = max(fx,0)
    fy = max(fy,0)
    fxy = max(fxy,0)

    # Compute final value and return
    xCy = fx - fy - fxy
    return xCy
    
def approx_logn_factorial(n):
    '''
    Use Stirling's approximation to compute log(n!).
    '''
    x = (n+0.5)*log(n) - n + 0.5*log(2*pi)
    return x
