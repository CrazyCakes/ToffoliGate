import numpy as np
import qutip as qt
#from tqdm import tqdm
N_qi=5 #number of information qubits
N=N_qi+2  #num of total qubits

n=0
w=np.array([1, 1, 1, 1])  # small indx first (alternating one first)
plateform="cir"
deco="on"

n_range=25
n_min=20
dis_n=np.arange(n_min, n_range, 4)

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

#J=2*np.pi*40*10**6      #2*np.pi*2*10**3        2*np.pi*40*10**6
J=J/np.max(np.abs(w))
omaga0=0#*np.e/2
#omaga=(-2*n+(N-1))*J-omaga0
omaga=-2*J-omaga0#(-2*n+np.sum(w))*J-omaga0
# Coupling strength and external magnetic field   #(0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111)

#T1=30*10**(-6)   #1  30*10**(-6)
#T2=30*10**(-6)      # 50  30*10**(-6)
#rt_gamma_amp=(1/T1)**0.5                  #  (2*np.pi*100*10**3)**0.5
#rt_gamma_phase=(abs(0.5*(1/T2-1/(2*T1))))**0.5                    #(1/T2)**0.5  #(2*np.pi*500*10**3)**0.5


# Construct the Pauli matrices


def Build_H_J(c1, c2, target):
    H_c1_target=qt.qeye(2)
    H_c2_target=qt.qeye(2)

    if c1==0:
        H_c1_target=qt.sigmaz()
    elif target==0:
        H_c1_target=qt.sigmaz()
    if c2==0:
        H_c2_target=qt.sigmaz()
    elif target==0:
        H_c2_target=qt.sigmaz()

    for i in range(1, N):
        if c1==i or target==i:
            H_c1_target=qt.tensor(qt.sigmaz(), H_c1_target)
        else:
            H_c1_target=qt.tensor(qt.qeye(2), H_c1_target)
        if c2==i or target==i:
            H_c2_target=qt.tensor(qt.sigmaz(), H_c2_target)
        else:
            H_c2_target=qt.tensor(qt.qeye(2), H_c2_target)
    return H_c1_target+H_c2_target

def Build_H_x(target):
    H_x=qt.qeye(2)
    if target==0:
        H_x=qt.sigmax()
    for i in range(1, N):
        if target==i:
            H_x=qt.tensor(qt.sigmax(), H_x)
        else:
            H_x=qt.tensor(qt.qeye(2), H_x)
    return H_x

def Build_H_y(target):
    H_y=qt.qeye(2)
    if target==0:
        H_y=qt.sigmay()
    for i in range(1, N):
        if target==i:
            H_y=qt.tensor(qt.sigmay(), H_y)
        else:
            H_y=qt.tensor(qt.qeye(2), H_y)
    return H_y



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
    
  


def osc_cos(t, args):
    return np.cos(omaga*t)
def osc_sin(t, args):
    return np.sin(omaga*t)

def osc_cos_pi(t, args):
    return np.cos(omaga*t+np.pi)
def osc_sin_pi(t, args):
    return np.sin(omaga*t+np.pi)


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

Ubases=np.empty((2**(2*N_qi)), dtype=object)
for i in range(2**N_qi): #|0><i|
    for j in range(2**N_qi): #prefactor
        if N_qi==3:
            U=basis_state(N_qi, 0)*basis_state(N_qi, i).dag()*prefactor_basis_3[j, 0]
        elif N_qi==4:
            U=basis_state(N_qi, 0)*basis_state(N_qi, i).dag()*prefactor_basis_4[j, 0]
        elif N_qi==5:
            U=basis_state(N_qi, 0)*basis_state(N_qi, i).dag()*prefactor_basis_5[j, 0]

        for k in range(1, 2**N_qi): # make each element
            if N_qi==3:
                U=U+prefactor_basis_3[j][k]*basis_state(N_qi, k)*basis_state(N_qi, (k+i)%(2**N_qi)).dag()
            elif N_qi==4:
                U=U+prefactor_basis_4[j][k]*basis_state(N_qi, k)*basis_state(N_qi, (k+i)%(2**N_qi)).dag()   
            elif N_qi==5:
                U=U+prefactor_basis_5[j][k]*basis_state(N_qi, k)*basis_state(N_qi, (k+i)%(2**N_qi)).dag()      
        Ubases[2**N_qi*i+j]=U#/((U.dag()*U).tr())**0.5

#######
#Toffoli
Toffoli=0
for i in range(2**(N_qi-1)):
    total_w=0
    temp_i=i
    for j in range(N_qi-1):
        if temp_i%2!=0:
            total_w=total_w+w[j]
        else:
            total_w=total_w-w[j]
        temp_i=temp_i//2
    if total_w==n:# 
        print(i)
        Toffoli=Toffoli+qt.tensor(basis_state(N_qi-1, i)*basis_state(N_qi-1, i).dag(), -1j*qt.sigmax())
    else:
        #print(i)
        #print(basis_state(N-1, i))
        #print(basis_state(N-1, i)*basis_state(N-1, i).dag())
        Toffoli=Toffoli+qt.tensor(basis_state(N_qi-1, i)*basis_state(N_qi-1, i).dag(), qt.qeye(2))



def X_gate(target):
    X=qt.qeye(2)
    if target==0:
        X=qt.sigmax()
    for i in range(1, N):
        if target==i:
            X=qt.tensor(qt.sigmax(), X)
        else:
            X=qt.tensor(qt.qeye(2), X)
    return X

def Toffoli_S(c1, c2, target, pi, Omaga, times):
    S=0
    if pi==0:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*Build_H_J(c1, c2, target), [Omaga*Build_H_x(target), osc_cos], [Omaga*Build_H_y(target), osc_sin]], tlist=times)
    else:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*Build_H_J(c1, c2, target), [Omaga*Build_H_x(target), osc_cos_pi], [Omaga*Build_H_y(target), osc_sin_pi]], tlist=times)
    if deco=="on":
        S=qt.liouvillian(S, [rt_gamma_amp*L_amp, rt_gamma_phase*L_phase])
    return S

def Two_Toffoli_S(c1_1, c1_2, target_1, c2_1, c2_2, target_2, pi, Omaga, times):
    S=0
    if pi==0:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*Build_H_J(c1_1, c1_2, target_1)+J/2*Build_H_J(c2_1, c2_2, target_2), [Omaga*Build_H_x(target_1), osc_cos], [Omaga*Build_H_y(target_1), osc_sin], [Omaga*Build_H_x(target_2), osc_cos], [Omaga*Build_H_y(target_2), osc_sin]], tlist=times)
    else:
        S=qt.QobjEvo([-omaga0/2*H_E0+J/2*Build_H_J(c1_1, c1_2, target_1)+J/2*Build_H_J(c2_1, c2_2, target_2), [Omaga*Build_H_x(target_1), osc_cos_pi], [Omaga*Build_H_y(target_1), osc_sin_pi], [Omaga*Build_H_x(target_2), osc_cos_pi], [Omaga*Build_H_y(target_2), osc_sin_pi]], tlist=times)
    if deco=="on":
        S=qt.liouvillian(S, [rt_gamma_amp*L_amp, rt_gamma_phase*L_phase])
    return S
    

ancilla=qt.tensor(qt.basis(2, 0)*qt.basis(2, 0).dag(), qt.basis(2, 0)*qt.basis(2, 0).dag())
Ubases_ancilla=np.empty((2**(2*N_qi)), dtype=object)
for i in range(2**(2*N_qi)): 
    Ubases_ancilla[i]=qt.tensor(ancilla, Ubases[i])

Toffoli_ancilla=qt.tensor(qt.tensor(qt.qeye(2), qt.qeye(2)),  Toffoli)

Ubases_ref_dag=np.empty((2**(2*N_qi)), dtype=object)
for i in range(2**(2*N_qi)): 
    Ubases_ref_dag[i]=Toffoli*Ubases[i].dag()*Toffoli.dag()

X_Gates=np.empty(N, dtype=object)
for i in range(N):
    X_Gates[i]=X_gate(i)

X_Gates_dag=np.empty(N, dtype=object)
for i in range(N):
    X_Gates_dag[i]=X_Gates[i].dag()

X_Pos=np.array([[3, 4], [2, 4], [2, 3], [1, 3], [1, 2], [1, 4]])
Toffoli_Gates=np.empty(7, dtype=object)
Toffoli_pi_Gates=np.empty(6, dtype=object)

file = open("./fidelity_cir_20.txt", 'w')
options=qt.Options(store_final_state=True)
for i in range(len(dis_n)):
    Omaga=J/(dis_n[i]) #8 is general
    T=np.pi/(2*Omaga)
    times=np.linspace(0, T, steps)
    Toffoli_Gates[0]=Two_Toffoli_S(1, 2, 5, 3, 4, 6, 0, Omaga, times)
    Toffoli_Gates[1]=Two_Toffoli_S(1, 3, 5, 2, 4, 6, 0, Omaga, times)
    Toffoli_Gates[2]=Two_Toffoli_S(1, 4, 5, 2, 3, 6, 0, Omaga, times)
    Toffoli_Gates[3]=Two_Toffoli_S(2, 4, 5, 1, 3, 6, 0, Omaga, times)
    Toffoli_Gates[4]=Two_Toffoli_S(3, 4, 5, 1, 2, 6, 0, Omaga, times)
    Toffoli_Gates[5]=Two_Toffoli_S(2, 3, 5, 1, 4, 6, 0, Omaga, times)

    Toffoli_Gates[6]=Toffoli_S(5, 6, 0, 0, Omaga, times)

    Toffoli_pi_Gates[0]=Two_Toffoli_S(1, 2, 5, 3, 4, 6, 1, Omaga, times)
    Toffoli_pi_Gates[1]=Two_Toffoli_S(1, 3, 5, 2, 4, 6, 1, Omaga, times)
    Toffoli_pi_Gates[2]=Two_Toffoli_S(1, 4, 5, 2, 3, 6, 1, Omaga, times)
    Toffoli_pi_Gates[3]=Two_Toffoli_S(2, 4, 5, 1, 3, 6, 1, Omaga, times)
    Toffoli_pi_Gates[4]=Two_Toffoli_S(3, 4, 5, 1, 2, 6, 1, Omaga, times)
    Toffoli_pi_Gates[5]=Two_Toffoli_S(2, 3, 5, 1, 4, 6, 1, Omaga, times)
    Fidelity=0
    print("start")
    for j in range(2**(2*N_qi)):   #tqdm(range(2**(2*N_qi))):
        #print("calculating")
        state=Ubases_ancilla[j]
        for k in range(6):
            state=X_Gates[X_Pos[k][0]]*X_Gates[X_Pos[k][1]]*state*X_Gates_dag[X_Pos[k][1]]*X_Gates_dag[X_Pos[k][0]]
            result=qt.mesolve(Toffoli_Gates[k], state, times, options=options)
            result=qt.mesolve(Toffoli_Gates[6], result.final_state, times, options=options)
            result=qt.mesolve(Toffoli_pi_Gates[k], result.final_state, times, options=options)
            state=X_Gates[X_Pos[k][0]]*X_Gates[X_Pos[k][1]]*result.final_state*X_Gates_dag[X_Pos[k][1]]*X_Gates_dag[X_Pos[k][0]]
        Fidelity=Fidelity+(Ubases_ref_dag[j]*state.ptrace([2, 3, 4, 5, 6])).tr() 
        #Fidelity=Fidelity+(Toffoli_ancilla*Ubases_ancilla[j].dag()*Toffoli_ancilla.dag()*state).tr()      
    print(str(J/Omaga)+" "+str((Fidelity+2**(2*N_qi))/(2**(2*N_qi)*(2**N_qi+1))))
    print(str(J/Omaga)+" "+str((Fidelity+2**(2*N_qi))/(2**(2*N_qi)*(2**N_qi+1))), file=file)
    file.flush()
