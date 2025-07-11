import numpy as np
import qutip as qt
N=4  #num of qubits
n1=1
n2=3
n3=2
w=np.array([1, 1, 1, 1])  # small indx first (alternating one first)
plateform="ion" #cir, ion
deco="on" #on, off
number_drive=2 #1, 2, 3
n_range=33
n_min=4
dis_n=np.arange(n_min, n_range, 4)
#dis_n=np.array([2, 4, 6, 8, 7, 13, np.e*7])
steps=100

J=0
T1=0
T2=0
rt_gamma_amp=0
rt_gamma_phase=0
if plateform=="cir":
    J=2*np.pi*40*10**6 
    T1=30*10**(-6)
    T2=30*10**(-6)
    rt_gamma_amp=(1/T1)**0.5
    rt_gamma_phase=(1/T2-1/(2*T1))**0.5
elif plateform=="ion":
    J=2*np.pi*2*10**3 
    T2=50
    rt_gamma_amp=0
    rt_gamma_phase=(1/T2)**0.5


J=J/np.max(np.abs(w))
omaga0=0#*np.e/2
#omaga=(-2*n+(N-1))*J-omaga0
omaga1=-n1*J-omaga0#(-2*n1+np.sum(np.abs(w)))*J-omaga0
omaga2=-n2*J-omaga0#(-2*n2+np.sum(np.abs(w)))*J-omaga0
omaga3=-n3*J-omaga0#(-2*n3+np.sum(np.abs(w)))*J-omaga0
# Coupling strength and external magnetic field   #(0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111)

#T1=1   #1   30*10**(-6)
#T2=50     # 50  30*10**(-6)
#rt_gamma_amp=0    #(1/T1)**0.5                  #  (2*np.pi*100*10**3)**0.5
#rt_gamma_phase=(0.5*(1/T2))**0.5         #(abs(0.5*(1/T2-1/(2*T1))))**0.5                    #(1/T2)**0.5  #(2*np.pi*500*10**3)**0.5


# Construct the Pauli matrices
H_J=0
# Create the Hamiltonian
for i in range(1):
    for j in range(i + 1, N):
        if i==0:
            pH=qt.sigmaz()*w[j-1]
        else:
            pH=qt.qeye(2)
        for k in range(1, N):
            if k==i or k==j:
                pH=qt.tensor(qt.sigmaz(), pH)
            else:
                pH=qt.tensor(qt.qeye(2), pH)   
        H_J=H_J+pH



H_E0=0
for i in range(N):
    if i==0:
        pH=qt.sigmaz()
    else:
        pH=qt.qeye(2)
    for j in range(1, N):
        if j==i:
            pH=qt.tensor(qt.sigmaz(), pH)
        else:
            pH=qt.tensor(qt.qeye(2), pH)
    H_E0=H_E0+pH
    
  


H_dx=qt.sigmax()
for i in range(N-1):
    H_dx=qt.tensor(qt.qeye(2), H_dx)
    #H_dx=qt.tensor(H_dx, qt.qeye(2))

def osc_cos1(t, args):
    return np.cos(omaga1*t)

def osc_cos2(t, args):
    return np.cos(omaga2*t)

def osc_cos3(t, args):
    return np.cos(omaga3*t)

H_dy=qt.sigmay()
for i in range(N-1):
    H_dy=qt.tensor(qt.qeye(2), H_dy)
    #H_dy=qt.tensor(H_dy, qt.qeye(2))


def osc_sin1(t, args):
    return np.sin(omaga1*t)

def osc_sin2(t, args):
    return np.sin(omaga2*t)

def osc_sin3(t, args):
    return np.sin(omaga3*t)


L_amp=0
for i in range(N):
    if i==0:
        pH=qt.basis(2, 0)*qt.basis(2, 1).dag()
    else:
        pH=qt.qeye(2)
    for j in range(1, N):
        if j==i:
            pH=qt.tensor(qt.basis(2, 0)*qt.basis(2, 1).dag(), pH)
        else:
            pH=qt.tensor(qt.qeye(2), pH)
    L_amp=L_amp+pH

L_phase=0
for i in range(N):
    if i==0:
        pH=qt.sigmaz()
    else:
        pH=qt.qeye(2)
    for j in range(1, N):
        if j==i:
            pH=qt.tensor(qt.sigmaz(), pH)
        else:
            pH=qt.tensor(qt.qeye(2), pH)
    L_phase=L_phase+pH


#making unitary bases

def basis_state(dim, label):    #2**3 deim
    state=0

    pstate=qt.basis(2, label%2)
    label=label//2
    for k in range(dim-1):
        pstate=qt.tensor(qt.basis(2, label%2), pstate)
        label=label//2
    state=state+pstate

    state=state.unit() 
    return state 

prefactor_basis_3= np.array([
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, -1, -1, -1, -1],
    [1, 1, -1, -1, 1, 1, -1, -1],
    [1, -1, 1, -1, 1, -1, 1, -1],
    [1, 1, -1, -1, -1, -1, 1, 1],
    [-1, 1, 1, -1, -1, 1, 1, -1],
    [1, -1, 1, -1, -1, 1, -1, 1],
    [-1, 1, 1, -1, 1, -1, -1, 1], 
    ])  

prefactor_basis_4= np.array([
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1],
    [1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1],
    [1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1],
    [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
    [1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1],
    [1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
    [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1],
    [-1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1],
    [-1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1],
    [1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1],
    [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1],
    [-1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1], 
    [-1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1]
    ])      

prefactor_basis_5= np.array([
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1],
    [1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1],
    [1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1],
    [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
    [1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1],
    [1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
    [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1],
    [-1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1],
    [-1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1],
    [1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1],
    [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1],
    [-1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1],
    [-1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1],
    [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1],
    [1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1],
    [1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1],
    [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1],
    [1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1],
    [1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1],
    [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1],
    [-1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1],
    [-1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1],
    [1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1],
    [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1],
    [-1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1],
    [-1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1]
    ])      

Ubases=np.empty((2**(2*N)), dtype=object)
for i in range(2**N): #|0><i|
    for j in range(2**N): #prefactor
        if N==3:
            U=basis_state(N, 0)*basis_state(N, i).dag()*prefactor_basis_3[j, 0]
        elif N==4:
            U=basis_state(N, 0)*basis_state(N, i).dag()*prefactor_basis_4[j, 0]
        elif N==5:
            U=basis_state(N, 0)*basis_state(N, i).dag()*prefactor_basis_5[j, 0]

        for k in range(1, 2**N): # make each element
            if N==3:
                U=U+prefactor_basis_3[j][k]*basis_state(N, k)*basis_state(N, (k+i)%(2**N)).dag()
            elif N==4:
                U=U+prefactor_basis_4[j][k]*basis_state(N, k)*basis_state(N, (k+i)%(2**N)).dag()   
            elif N==5:
                U=U+prefactor_basis_5[j][k]*basis_state(N, k)*basis_state(N, (k+i)%(2**N)).dag()      
        Ubases[2**N*i+j]=U#/((U.dag()*U).tr())**0.5

#######
#Toffoli
Toffoli=0
if number_drive==1:
    for i in range(2**(N-1)):
        total_w=0
        temp_i=i
        for j in range(N-1):
            if temp_i%2!=0:
                total_w=total_w+w[j]
            else:
                total_w=total_w-w[j]
            temp_i=temp_i//2
        if total_w==n1:# 
            print(i)
            Toffoli=Toffoli+qt.tensor(basis_state(N-1, i)*basis_state(N-1, i).dag(), -1j*qt.sigmax())
        else:
            Toffoli=Toffoli+qt.tensor(basis_state(N-1, i)*basis_state(N-1, i).dag(), qt.qeye(2))
elif number_drive==2:
    for i in range(2**(N-1)):
        total_w=0
        temp_i=i
        for j in range(N-1):
            if temp_i%2!=0:
                total_w=total_w+w[j]
            else:
                total_w=total_w-w[j]
            temp_i=temp_i//2
        if total_w==n1 or total_w==n2:# 
            print(i)
            Toffoli=Toffoli+qt.tensor(basis_state(N-1, i)*basis_state(N-1, i).dag(), -1j*qt.sigmax())
        else:
            Toffoli=Toffoli+qt.tensor(basis_state(N-1, i)*basis_state(N-1, i).dag(), qt.qeye(2))
elif number_drive==3:
    for i in range(2**(N-1)):
        total_w=0
        temp_i=i
        for j in range(N-1):
            if temp_i%2!=0:
                total_w=total_w+w[j]
            else:
                total_w=total_w-w[j]
            temp_i=temp_i//2
        if total_w==n1 or total_w==n2 or total_w==n3:# 
            print(i)
            Toffoli=Toffoli+qt.tensor(basis_state(N-1, i)*basis_state(N-1, i).dag(), -1j*qt.sigmax())
        else:
            Toffoli=Toffoli+qt.tensor(basis_state(N-1, i)*basis_state(N-1, i).dag(), qt.qeye(2))


options=qt.Options(store_final_state=True)
for i in range(len(dis_n)):
    Omaga=J/(dis_n[i]) #8 is general
    T=np.pi/(2*Omaga)
    times=np.linspace(0, T, steps)
    S=0
    if number_drive==1:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos1], [Omaga*H_dy, osc_sin1]], tlist=times)
    elif number_drive==2:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos1], [Omaga*H_dx, osc_cos2], [Omaga*H_dy, osc_sin1], [Omaga*H_dy, osc_sin2]], tlist=times)
    elif number_drive==3:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos1], [Omaga*H_dx, osc_cos2], [Omaga*H_dx, osc_cos3], [Omaga*H_dy, osc_sin1], [Omaga*H_dy, osc_sin2], [Omaga*H_dy, osc_sin3]], tlist=times)
    if deco=="on":
        S=qt.liouvillian(S, [rt_gamma_amp*L_amp, rt_gamma_phase*L_phase])
    #H_1=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos1], [Omaga*H_dy, osc_sin1]], tlist=times)
    #H_2=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos1], [Omaga*H_dx, osc_cos2], [Omaga*H_dy, osc_sin1], [Omaga*H_dy, osc_sin2]], tlist=times)
    #H_3=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos1], [Omaga*H_dx, osc_cos2], [Omaga*H_dx, osc_cos3], [Omaga*H_dy, osc_sin1], [Omaga*H_dy, osc_sin2], [Omaga*H_dy, osc_sin3]], tlist=times)
    #H_mul=qt.QobjEvo([-omaga0/2*H_E0+J/2*H_J, [Omaga*H_dx, osc_cos12], [Omaga*H_dy, osc_sin12]], tlist=times)
    #S_1=qt.liouvillian(H_1, [rt_gamma_amp*L_amp, rt_gamma_phase*L_phase])
    #S_2=qt.liouvillian(H_2, [rt_gamma_amp*L_amp, rt_gamma_phase*L_phase])
    #S_3=qt.liouvillian(H_3, [rt_gamma_amp*L_amp, rt_gamma_phase*L_phase])
    #H_rJ=qt.QobjEvo([-omaga0/2*H_E0-J/2*H_J], tlist=times)
    Fidelity=0
    for j in range(2**(2*N)):
        #print(Ubases[j])
        #print((Ubases[j]*Ubases[j].dag()).tr())
        #result=qt.sesolve(H, qt.tensor(qt.basis(2, 1), qt.basis(2, 1), qt.basis(2, 1)), times, [], options=options)
        #print(result.final_state)
        result=qt.mesolve(S, Ubases[j], times, options=options)
        #result=qt.mesolve(H_rJ, result.final_state, times, options=options)
        #result=qt.mesolve(H_mul, result.final_state, times, options=options)
        #result=qt.mesolve(H, Ubases[j], times, options=options)
        #result=qt.mesolve(H_rJ, result.final_state, times, options=options)
        #Fidelity=Fidelity+(qt.toffoli()*Ubases[j].dag()*qt.toffoli().dag()*result.final_state).tr()
        #Fidelity=Fidelity+(I*Ubases[j].dag()*I.dag()*result.final_state).tr()
        Fidelity=Fidelity+(Toffoli*Ubases[j].dag()*Toffoli.dag()*result.final_state).tr()
        
        #Fidelity=Fidelity+(Toffoli*Ubases[j].dag()*Toffoli.dag()*Toffoli*Ubases[j]).tr()
        #Fidelity=Fidelity+(I*Ubases[j].dag()*I.dag()*I*Ubases[j]).tr()
        #print(result.final_state)
    #print(Fidelity)
    print(str(((Fidelity+2**(2*N))/(2**(2*N)*(2**N+1))).real))
    #print(str(J/Omaga)+" "+str(((Fidelity+2**(2*N))/(2**(2*N)*(2**N+1))).real))
