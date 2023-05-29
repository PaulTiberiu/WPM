import numpy as np
import matplotlib.pyplot as plt

#1

N = 1000
c0 = 343
rho0 = 1
l = 0.35 #m
x = np.linspace(0,l,1000)
#print(x)
j = 1
h = np.zeros(1000)
S = np.zeros(1000)
k = 0
S0=0.06 #m
Sf=0.18 #m
f = np.linspace(1, 10000, 1000)
#print(f)

for j in range(x.size):
    h[j] = x[j]-x[j-1]
    
h[0] = 0.0004004

#0.45 pour l=0.4
#avant 3,14
beta = 0.005
for k in range (S.size):
    S[k] = S0*np.exp(beta*x[k])
"""
x0 = 0.14
for k in range(S.size):
    S[k] = S0*(x[k]**2/x0**2)
"""
fc = (beta*c0)/(4*np.pi)
print("freq coupure", fc)
#print(S)beta

#plt.plot(x, S)
#plt.show()

def Tn(f0,rho0,c0,S,h):
    Z0 = rho0*c0
    k0 = 2*np.pi*f0/c0
    tn = np.zeros((2,2),dtype=complex) #matrice 2,2 complexe init a 0
    tn[0,0] = np.cos(k0*h)
    tn[0,1] = -1j*Z0/S*np.sin(k0*h)
    tn[1,0] = -1j*S/Z0*np.sin(k0*h)
    tn[1,1] = np.cos(k0*h)
    return tn

D0 = np.zeros((2,2))
Z0 = rho0*c0
D0[0,0] = 1
D0[0,1] = 1
D0[1,0] = S0/Z0
D0[1,1] = -S0/Z0
D0_inv = np.linalg.inv(D0)
#M = np.zeros((2,2),dtype=complex)
M = D0_inv

for i in range(N):
    M = np.dot(M,Tn(f[0],rho0,c0,S[i],h[i]))

Df = np.zeros((2,2))
Df[0,0] = 1
Df[0,1] = 1
Df[1,0] = Sf/Z0
Df[1,1] = -Sf/Z0

M = np.dot(M,Df)

#4
t = 1/M[0,0]
r = M[1,0]*t
print("le coefficient t est: ",t)
print("le coefficient r est: ",r)



#print(np.abs(r))
R = (np.abs(r))**2
k0 = 2*np.pi*f[0]/c0
T = (np.abs(t))**2
#print("R : ",R)
#print("T: ",T)
I_i= 10**(-7)
L_i = 64
I_t = T*I_i
#print("I_t: ",I_t)
L_t = 10*np.log10(I_t/10**(-12))
#print("L_t: ",L_t)


#5
def Mn(f):
    D0 = np.zeros((2,2))
    Z0 = rho0*c0
    D0[0,0] = 1
    D0[0,1] = 1
    D0[1,0] = S0/Z0
    D0[1,1] = -S0/Z0
    D0_inv = np.linalg.inv(D0)
    M = D0_inv

    for i in range(N):
        M = np.dot(M,Tn(f,rho0,c0,S[i],h[i]))

    Df = np.zeros((2,2))
    Df[0,0] = 1
    Df[0,1] = 1
    Df[1,0] = Sf/Z0
    Df[1,1] = -Sf/Z0

    M = np.dot(M,Df)
    return M

t = np.zeros(1000,dtype=complex) #tableau de coef transmission
for i in range(1000):
    M = Mn(f[i])
    t[i] = 1/M[0,0]
#print(t)
#plt.plot(f,t)
#plt.show()
#print("M[1,1] = ",M[0,0])
#print("M[2,2] = ",M[1,1])
#print("M[1,2] = ",M[0,1])
#print("M[2,1] = ",M[1,0])

#6
#tracer facteur de perte en fonction de la frequence
facteur_perte = 10*np.log10(1/(np.abs(t)**2))
plt.plot(f,facteur_perte)
plt.xlabel('$f$ (Hz)')
plt.ylabel('$L_T$ (dB)')
plt.xlim(0,4000)
plt.show()