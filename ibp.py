# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as st
import scipy.misc as misc
from scipy.special import gammaln
from collections import Counter

def buffet(alpha,N):
    Z=np.ones(shape=(1,max(1,st.poisson(alpha).rvs())))
    for i in xrange(2,N+1):
        old=np.array([[st.bernoulli(np.sum(Z[:,k])/float(i)).rvs() for k in xrange(len(Z[0]))]])
        new=np.ones(shape=(1,st.poisson(alpha/float(i)).rvs()))
        new_cols=np.vstack((np.zeros(shape=(i-1,len(new[0]))),new))
        Z=np.hstack((np.vstack((Z,old)),new_cols))
    return Z

def likelihood(alpha,Z):
    N,K=Z.shape
    bin_counts=Counter(map(tuple,Z.T))
    all_bin=[tuple([float(b) for b in "{0:b}".format(i+1)]) for i in xrange(0,2**Z.shape[0]-1)]
    column_reps=np.sum([misc.factorial(bin_counts[b]) if b in bin_counts.keys() else 1 for b in all_bin])
    
    dispersion=(alpha**K)/column_reps
    harm=np.exp(-alpha*np.sum([1./n for n in xrange(1,N+1)]))
    count=np.prod([(misc.factorial(N-m)*misc.factorial(m-1))/misc.factorial(N) for m in Z.sum(axis=0)])
    return dispersion*harm*count
    
def noisy_or(eps,lamb,Y,Z,i,t,x=1):
    l=1-((1-lamb)**np.dot(Z[i,:],Y[:,t])*(1-eps))
    return (l**x)-(l*(1-x))
    
def cond_noisy_orZ(eps,lamb,Y,Z,i,t,k,z,x):
    propZ=Z.copy()
    propZ[i][k]=z
    return noisy_or(eps,lamb,Y,propZ,i,t,x)
    
def cond_noisy_orY(eps,lamb,Y,Z,i,t,k,y,x):
    propY=Y.copy()
    propY[k][t]=y
    return noisy_or(eps,lamb,propY,Z,i,t,x)
    
def marginal_like(eps,lamb,X,Y,Z):
    N,T=X.shape
    return np.prod([[noisy_or(eps,lamb,Y,Z,i,t,X[i][t]) for i in xrange(N)] for t in xrange(T)])
    
def log_like(alpha,p,eps,lamb,X,Y,Z):
    prior=np.log(likelihood(alpha,Z))
    like=np.log(marginal_like(eps,lamb,X,Y,Z))
    return prior+like
    
def cond_new(p,eps,lamb,X,Y,Z,i,limit=10):
    l=np.ones(limit)
    T=X.shape[1]
    for t in xrange(T):
        eta=1-(1-eps)**np.dot(Z[i],Y[:,t])
        l*=np.array([1-(1-eps)*eta*(1-lamb*p)**Knew for Knew in xrange(1,limit+1)])
    return l
    
def sample_oldZ(eps,lamb,X,Y,Z):
    N,K=Z.shape
    m=Z.sum(axis=0)
    for i in xrange(N):
        for k in xrange(K):
            mk=m[k]-Z[i][k]
            if mk>0:
                prior=mk/float(N)
                like0,like1=[np.prod(v) for v in zip(*np.array([[cond_noisy_orZ(eps,lamb,Y,Z,i,t,k,z,x=X[i][t]) for z in (0,1)] for t in xrange(T)]))]
                post=(prior*like1)/float(prior*like1+(1-prior)*like0)
                if np.random.uniform(0,1)<post:
                    Z[i][k]=1
                else:
                    Z[i][k]=0
    return Z
    
def sample_newZ(alpha,p,eps,lamb,X,Y,Z,limit=10):
    N,K=Z.shape
    m=Z.sum(axis=0)
    prior=[st.poisson(alpha/float(N)).pmf(ki) for ki in xrange(1,limit+1)]
    all_new=[]
    totalK=0
    to_zero=[]
    for i in xrange(N):
        for k in xrange(K):
            mk=m[k]-Z[i][k]
            if mk==0:
                if i not in [tz[0] for tz in to_zero]:
                    like=cond_new(p,eps,lamb,X,Y,Z,i,limit=10)
                    post=(prior*like)/np.sum(prior*like)
                    Knew=np.random.choice(range(1,limit+1),p=post)
                    totalK+=Knew
                    new_columns=np.zeros(shape=(N,Knew))
                    new_columns[i]=np.ones(Knew)
                    all_new.append(new_columns)
                to_zero.append((i,k))
    for i,k in to_zero:
        Z[i][k]=0
    if len(all_new)>0:
        all_new=np.hstack(all_new)
        Z=np.hstack((Z,all_new))
    return Z
    
def sample_Y_post(eps,lamb,p,X,Y,Z,k,t):
    like0,like1=[np.prod(v) for v in zip(*np.array([[cond_noisy_orY(eps,lamb,Y,Z,i,t,k,z,x=X[i][t]) for z in (0,1)] for i in xrange(N)]))]
    post=(p*like1)/float(p*like1+(1-p)*like0)
    if np.random.uniform(0,1)<post:
        return 1
    else:
        return 0
    
    
def sample_Y(eps,lamb,p,X,Y,Z):
    N,newK=Z.shape
    K,T=Y.shape
    Y=np.vstack((Y,np.random.choice([0,1],size=(newK-K,T))))
    for k in xrange(K+1,newK):
        for t in xrange(T):
            Y[k][t]=sample_Y_post(eps,lamb,p,X,Y,Z,k,t)
    for k in xrange(newK):
        for t in xrange(T):
            Y[k][t]=sample_Y_post(eps,lamb,p,X,Y,Z,k,t)
    return Y
    
def sample_post_alpha(Z,maxtries=100):
    proposal_dist=st.gamma(1,1)
    theta=np.linspace(*proposal_dist.interval(.95))
    f=proposal_dist.pdf(theta)
    g=likelihood(theta,Z)
    ratio=f/g
    M=theta[list(ratio).index(max(ratio))]
    tries=0
    while tries<maxtries:
        tries+=1
        y=proposal_dist.rvs()
        u=np.random.uniform(0,1)
        fy=proposal_dist.pdf(y)
        gy=likelihood(y,Z)
        print fy,gy,M,fy/(M*gy),u
        if u<fy/(M*gy):
            return y
            
def logPXYZ(alpha,p,eps,lamb,X,Y,Z):
    N,K=Z.shape
    T=X.shape[1]
    Hn=np.sum([1./i for i in xrange(1,N+1)])
    mk = np.sum(Z,axis=0)
    lPZ = K*np.log(alpha) -alpha*Hn + np.sum(gammaln(mk)+gammaln(N-mk+1)-gammaln(N+1))
    lPX = np.sum([[noisy_or(eps,lamb,Y,Z,i,t,x=X[i][t]) for i in xrange(N)] for t in xrange(T)])
    #print lPZ, lPX
    return lPZ+lPX

def eps_metropolis(alpha,p,eps,lamb,prop_variance,accepted,iteration,X,Y,Z):
    prop_eps=eps+np.random.randn()*prop_variance
    if prop_eps>0 and prop_eps<1:
        lold=logPXYZ(alpha,p,eps,lamb,X,Y,Z)
        lnew=logPXYZ(alpha,p,prop_eps,lamb,X,Y,Z)
        ratio=lnew-lold
        u=np.log(np.random.uniform(0,1))
        if u<min(0,ratio):
            eps=prop_eps
            accepted+=1
    acceptance_ratio=accepted/float(iteration)
    if iteration>20 and iteration%5==0:
        if acceptance_ratio<.2:
            prop_variance*=.9
        elif acceptance_ratio>.3:
            prop_variance*=1.1
        if prop_variance>.5:
            prop_variance=.5
    return eps,prop_variance,accepted
    
def lamb_metropolis(alpha,p,eps,lamb,prop_variance,accepted,iteration,X,Y,Z):
    prop_lamb=lamb+np.random.randn()*prop_variance
    if prop_lamb>0 and prop_lamb<1:
        lold=logPXYZ(alpha,p,eps,lamb,X,Y,Z)
        lnew=logPXYZ(alpha,p,eps,prop_lamb,X,Y,Z)
        ratio=lnew-lold
        u=np.log(np.random.uniform(0,1))
        if u<min(0,ratio):
            lamb=prop_lamb
            accepted+=1
    acceptance_ratio=accepted/float(iteration)
    if iteration>20 and iteration%5==0:
        if acceptance_ratio<.2:
            prop_variance*=.9
        elif acceptance_ratio>.3:
            prop_variance*=1.1
        if prop_variance>.5:
            prop_variance=.5
    return lamb,prop_variance,accepted
                         
def gibbs(alpha,p,eps,lamb,X,iterations=1,limit=10):
    eps_prop_variance=.05
    eps_accepted=0
    lamb_prop_variance=.05
    lamb_accepted=0
    
    N,T=X.shape
    Hn=np.sum([1./i for i in xrange(1,N+1)])
    Z=buffet(alpha,N)
    K=Z.shape[1]
    Y=np.random.choice([0,1],p=(1-p,p),size=(K,T))
    mode=0
    Zmode=None
    Ymode=None
    for iteration in xrange(iterations):
        Z=sample_oldZ(eps,lamb,X,Y,Z)
        Z=sample_newZ(alpha,p,eps,lamb,X,Y,Z,limit=10)
        Y=sample_Y(eps,lamb,p,X,Y,Z)
        nonzero=np.nonzero(Z.sum(axis=0))[0]
        Z=Z[:,nonzero]
        Y=Y[nonzero]
        K=Z.shape[1]
        l=marginal_like(eps,lamb,X,Y,Z)
        if l>mode:
            Zmode=Z
            Ymode=Y
        
        eps,eps_prop_variance,eps_accepted=eps_metropolis(alpha,p,eps,lamb,eps_prop_variance,eps_accepted,iteration+1,X,Y,Z)
        lamb,lamb_prop_variances,lamb_accepted=lamb_metropolis(alpha,p,eps,lamb,lamb_prop_variance,lamb_accepted,iteration+1,X,Y,Z)
        alpha=np.random.gamma(1+K,1/(1+Hn))
        ones=np.sum(Y)
        p=np.random.beta(ones+1,Y.size-ones+1)
    return Ymode,Zmode,alpha,p,eps,lamb

        
alpha=1
p=.5
eps=.01
lamb=.5
    
trueZ=np.array([[0,1,0],
                [1,0,0],
                [0,0,1]])
N,K=trueZ.shape
T=100
trueY=np.array([[1 if np.random.uniform(0,1)<p else 0 for t in xrange(T)] for k in xrange(K)])
X=np.array([[1 if np.random.uniform(0,1)<noisy_or(eps,lamb,trueY,trueZ,i,t) else 0 for t in xrange(T)] for i in xrange(N)])

Y,Z,alpha,p,eps,lamb=gibbs(alpha,p,eps,lamb,X,100)
