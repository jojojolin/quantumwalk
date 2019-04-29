#!/usr/bin/python
import numpy as np
import math
import matplotlib.pyplot as plt



def quantumwalk(steps, cstate, pos):
    t = 0 #counter
    # Initialise matrices
    k = np.zeros((steps, 2**steps),dtype=np.complex_)#complex projective space
    z = np.zeros((steps, 2**steps)) #position space
    y = np.zeros((steps, 2**steps)) #coin space
    
    # Hadamard coin (asymm)
    h = [[1/math.sqrt(2),1/math.sqrt(2)],[1/math.sqrt(2),-1/math.sqrt(2)]]
    

    # Grover's coin (symm)
    g = [[1/math.sqrt(2),1j/math.sqrt(2)],[1j/math.sqrt(2),1/math.sqrt(2)]]
    
    # Uncomment this part if symmetrical coin is preferred
    """
    h = g
    """
    for r in xrange(steps):
        for s in xrange(2**(r+1)):
            # Setting the coin states for all position states of the particle/walker
            if (s+1) % 2==0: #if the position is even
                y[r][s] = 1
            else:
                y[r][s] = 0
            # Setting the complex conjugates of the initial position states
            if r==0: #if it is the first iteration
                if cstate==0: #if the initial state of coin is 0
                    if s==0: #The very first element of the matrix
                        k[r][s] = h[0][0] # it will have value a
                    else:
                        k[r][s] = h[0][1]
                elif s==0: #if the initial state of coin is 1
                    k[r][s] = h[1][0] #The very first element of the matrix will have value c
                else:
                    k[r][s] = h[1][1] #the rest of the elt in the row will have value d
            else: #If it is not the first iteration 
                t=t+1
            if t==1: #coin flip transition
                k[r][s] = h[0][0]*k[r-1][int(round(s/2))]
            elif t==2:
                k[r][s] = h[0][1]*k[r-1][int(round(s/2))]
            elif t==3:
                k[r][s] = h[1][0]*k[r-1][int(round(s/2))]
            elif t==4:
                k[r][s] = h[1][1]*k[r-1][int(round(s/2))]
            	
            if t==4:# Reset
                t=0
            if r==0: #first iteration
                if y[r][s]==0: #if the coin state of initial state is 0
                    z[r][s]=pos-1 #shift to the left
                else:
                    z[r][s]=pos+1 #shift to the right
            elif y[r][s]==0:
                z[r][s]=z[r-1][int(round(s/2))]-1
            else:
                z[r][s]=z[r-1][int(round(s/2))]+1
    """ -----DEBUG-----
    print('y is')
    print(y)
    print('k is')
    print(k)
    print('z is')
    print(z)
    -----"""

    [m, _bin] = histc(z[steps-1][:], np.unique(z[steps-1][:]))
    
    
    
    """-----DEBUG-----
    print('z is :')
    print(z[steps-1][:])
    print('m is')
    print(m)
    print('_bin is')
    print(_bin)
    -----"""
    
    multiplez= indices(m, lambda x: x > 1) #m[np.where(m>1)]
    
    """-----DEBUG-----
    print('multiplez is')
    print(multiplez)
    -----"""
    
    indexz = indices(ismember(_bin,multiplez+1), lambda x: x != 0)
    
    nz = len(indexz)
    
    """-----DEBUG-----
    print("indexz")
    print(indexz)
    print("nz is %s" %nz)
    -----"""
    
    if indexz.ndim > 1:
        bz = len(indexz[0])
    else:
        bz = 1
    
    """-----DEBUG-----
    print("bz is %s " %bz)
    -----"""
    
    mz = len(multiplez)
    if multiplez.ndim > 1:
        bz = len(multiplez[0])
    else:
        bz = 1

    rz = 2
    zz = np.zeros((steps,2**steps))
    zz[steps-1][0] = z[steps-1][0]
    zz[steps-1][1] = z[steps-1][2**steps-1]
    yy = np.zeros((steps,2**steps))
    yy[steps-1][0]=y[steps-1][0]
    yy[steps-1][1]=y[steps-1][2**steps-1]
    kk = np.zeros((steps,2**steps),dtype=np.complex_)
    kk[steps-1][0]=k[steps-1][0]
    kk[steps-1][1]=k[steps-1][2**steps-1]
    rz1=rz
    
    for r in range(1,2**steps-1):
        zz[steps-1][rz1] =z[steps-1][r]
        yy[steps-1][rz1] =y[steps-1][r]
        kk[steps-1][rz1] =k[steps-1][r]
        rz1 += 1
    rzz = rz
    pk = 0
    qk = 0
    
    kkz = []
    
    """-----DEBUG-----
    print("kk")
    print(kk)
    print('mz is %s' %mz)
    -----"""
    
    for i in xrange(mz):
        indexz1 = indices(ismember(_bin, multiplez[i]+1), lambda x: x != 0)#find(ismember(_bin, multiplez[i]))
        
        """-----DEBUG-----
        print("indexz1")
        print(indexz1)
        -----"""
        
        kkz.append(z[steps-1][indexz1[i]])
    
    """-----DEBUG-----
    print("kkz")
    print(kkz) # is the pos value that have multiple coefficients
    -----"""
    
    for i in xrange(mz):
        for j in range(rz, 2**steps):
            if yy[steps-1][j]==0:
                if zz[steps-1][j]==kkz[i]:
                    pk+=kk[steps-1][j]
            elif yy[steps-1][j]==1:
                if zz[steps-1][j]==kkz[i]:
                    qk+=kk[steps-1][j]
        zz[steps-1][rzz]=kkz[i]
        yy[steps-1][rzz]=0
        kk[steps-1][rzz]=pk
        rzz = rzz+1
        zz[steps-1][rzz]=kkz[i]
        yy[steps-1][rzz]=1
        kk[steps-1][rzz]=qk
        pk=0
        qk=0
        rzz=rzz+1
    
    """-----DEBUG-----
    print('zz is')
    print(zz)
    -----"""
    
    print('States after %s iterations of Coin-flip and Shift transformation:' % steps)
    p = 0
    k2 = np.zeros((steps,steps*2),dtype=np.complex_)
    
    """-----DEBUG-----
    print('rzz is %s' %rzz)
    -----"""
    
    for s in xrange(rzz):
        print('%1.2f + %1.2fj |%1.0f,%1.0f>'%(kk[steps-1][s].real,kk[steps-1][s].imag,zz[steps-1][s],yy[steps-1][s]))
        k2[steps-1][s]=(kk[steps-1][s])**2 #calculate probability
        p=p+k2[steps-1][s]

    #Presenting the outcome 
    print('\nTotal Probability = %1.4f + %1.4fj \n' % (p.real, p.imag))
    print('Normalisation check:')
    
    if round(p.real)==1 and round(p.imag)==0:
        print("true\n")
    else:
        print("false, p=%s\n"% abs(p))
    kkz=rz
    nrzz=(rzz-2)/2+2
    
    """-----DEBUG-----
    print('nrzz is %s' %nrzz)
    -----"""
    pzz = np.zeros((steps,nrzz))
    pzz[steps-1][0]=zz[steps-1][0]
    pzz[steps-1][nrzz-1]=zz[steps-1][1]
    pk2 = np.zeros((steps,nrzz))
    pk2[steps-1][0]=abs(kk[steps-1][0])**2	
    pk2[steps-1][nrzz-1]=abs(kk[steps-1][1])**2

    for s in range(1,nrzz-1):
        kkz = kkz+1
        pzz[steps-1][s]=zz[steps-1][kkz-1]
        pk2[steps-1][s]=abs(kk[steps-1][kkz-1])**2 +abs(kk[steps-1][kkz])**2
        kkz=kkz+1
  

    print('+=================================+')
    print('|Position Probability Distribution|')
    print('+---------------------------------+')
    
    for s in xrange(nrzz):
        print(' %1.0f %1.8f'%(pzz[steps-1][s],pk2[steps-1][s]))
    
    plt.plot(pzz[steps-1][0:nrzz],pk2[steps-1][0:nrzz])
    plt.ylabel('Probability')
    plt.xlabel('Position')
    plt.show()


def ismember(A, B):
    """Returns an array of 1s and 0s with the size of a, where a and b are both type np.arrays."""
    return np.in1d(A,B)

def histc(x, bins):
    """A python version of the matlab histc"""
    map_to_bins = np.digitize(x, bins)
    r = np.zeros(bins.shape)
    for i in map_to_bins:
        r[i-1] +=1
    return [r, map_to_bins]
    
def indices(a, func):
    """Returns indices of a where the element meets the condition specified by the func parameter"""
    return np.array([i for (i, val) in enumerate(a) if func(val)])

def main():
    print('Number of steps:')
    steps = input()
    print('Initial coin state: ')
    cstate = input()
    print('Initial position: ')
    pos = input()
    print('Localise the particle to state: |%1.0f,%1.0f>\n' % (pos,cstate))
    # Catch for edge case
    if steps == 0:
        print('+=================================+')
        print('|Position Probability Distribution|')
        print('+---------------------------------+')
        print(' %1.0f %1.8f'%(pos,1))
    else:
        quantumwalk(steps,cstate,pos)

if __name__ == "__main__":
    main()