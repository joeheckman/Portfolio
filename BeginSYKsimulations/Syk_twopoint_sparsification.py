#current time: 1.234 seconds for N=10. Sparsification is slower at low N (my computer) and faster at larger N. That's decent news. 
#Last test showed that generating the Hamiltonian for larger N (22 on 1.13.25) is very doable. 
#However, two-point calculation for this Hamiltonian is still much too slow. 
#By the way, Mathematica is 1000 times faster and simpler for sparse matrix multiplication. 
import numpy as np
import scipy.linalg as linalg
from scipy.sparse import csr_matrix,vstack,kron, save_npz
from scipy.sparse.linalg import expm as sparse_expm
import matplotlib.pyplot as plt
import time
start_time = time.time()
N_maj = 22

##find dtype=complex is the homework problem
op_size = 2**(N_maj/2)
print([f'size of operator is {op_size} x {op_size}'])
chiEven_loop = np.zeros((int(N_maj/2),int(op_size),int(op_size)),dtype=complex)
chiOdd_loop = np.zeros((int(N_maj/2),int(op_size),int(op_size)),dtype=complex)

##pauli matrices
paulix = np.array([[0,1],[1,0]])
pauliy = np.array([[0,-1*1j], [1*1j,0]])
pauliz = np.array([[1,0],[0,-1]])
identity = np.array([[1,0],[0,1]])

#there will always be (N_maj/2)-1 pauliz in front of the unique pauli
#then (N_maj/2)-1 identity after. 
left_petal = [pauliz for _ in range(int(N_maj/2)-1)]
right_petal = [identity for _ in range(int(N_maj/2)-1)]
even_pistil = [paulix for _ in range(1)] 
odd_pistil = [pauliy for _ in range(1)]

left_petal_sparse = [csr_matrix(pauliz) for _ in range(int(N_maj/2)-1)]
right_petal_sparse  = [csr_matrix(identity) for _ in range(int(N_maj/2)-1)]


chi_even_builder = [*left_petal, *even_pistil, *right_petal] ## a flower
chi_odd_builder = [*left_petal, *odd_pistil, *right_petal]   ## another one, different species, yellower



mat_pos = -1
for i in range(int(N_maj/2)-1,-1,-1):
    mat_pos += 1
    print(i)    
    count = 0
    temp = np.kron(chi_even_builder[i],chi_even_builder[i + 1]) ##nucleus of the Majorana operator
    temp_odd = np.kron(chi_odd_builder[i],chi_odd_builder[i + 1])
    
    while True:
        count += 1
        temp = np.kron(temp,chi_even_builder[i+1+count]); ##this satisfies N_maj = 6 and above
        temp_odd = np.kron(temp_odd,chi_odd_builder[i+1+count])
        print(i+1+count)
        
        if len(temp) == op_size:
            print(op_size)
            chiEven_loop[mat_pos,:,:] = temp;
            chiOdd_loop[mat_pos,:,:] = temp_odd;
            break

        
Hamiltonian = np.zeros((int(op_size),int(op_size)),dtype=complex)
sparse_Hamiltonian = csr_matrix((int(op_size),int(op_size)), dtype=complex)

mean = 0
std_dev = np.sqrt(((6))/(16*N_maj**3))

##make arrays of fermions 
Fermions = np.vstack((chiEven_loop,chiOdd_loop))
##now make it sparse 
sparse_Fermions = []

for n in range(0,N_maj):    
    sparse_Fermions.append(csr_matrix(Fermions[n,:,:]))
    

##classic 4 interaction case. can use a gigantic computer to simulate a huge number of interaction terms with very large operators. 
for i in range(0,int(N_maj)):
    for j in range(0,int(N_maj)):
        for k in range(0,int(N_maj)):
            for l in range(0,int(N_maj)):
                print(f'Chi{i}*Chi{j}*Chi{k}*Chi{l}')
                randomNumber = np.random.normal(mean,std_dev)
                
                sparse_Hamiltonian += randomNumber*(sparse_Fermions[i] @ sparse_Fermions[j] @ sparse_Fermions[k] @ sparse_Fermions[l])
endtime = time.time()

##this is for 4 interaction term, (i^(interaction terms/2))
sparse_Hamiltonian = -1*sparse_Hamiltonian
# diving into the 2point here 
beta = 1
range_values = np.arange(0,beta,.01) ##Periodicity at beta right?
array_length = len(range_values)
new_array2 = csr_matrix((int(N_maj),array_length))
new_array2Test = csr_matrix((int(N_maj),array_length))

imaginary_time_evo = sparse_expm(-beta*sparse_Hamiltonian) 
#for i in range(0,int(N_maj)): 
i = 1;
count = -1
for tau in np.arange(0,beta,.01):
    print(tau)
    count += 1
    
    ##1.13.25 gets totally stuck here. 
    
    
    #time_evo_operator is its own thing to keep all of the time evolution parts attached to a Majorana operator
    time_evo_operator = csr_matrix(sparse_expm((-beta+tau)*sparse_Hamiltonian) @ Fermions[i] @ sparse_expm(-1*sparse_Hamiltonian*tau))  ##they win again, that being the youtube video Leo and Evan keep telling me to watch and GPT
    #point = np.trace(time_evo_operator @ Fermions[i,:,:]) / np.trace(linalg.expm((-beta+tau)*Hamiltonian) @ linalg.expm(-1*Hamiltonian*tau))
    point = (time_evo_operator @ sparse_Fermions[i]) / (sparse_expm((-beta+tau)*sparse_Hamiltonian) @ sparse_expm(-1*sparse_Hamiltonian))    
    point_trace = point.trace()
    
    new_array2[i,count]=point_trace

twoPoint = np.mean(new_array2,axis=0)
totaltime = endtime-start_time
print(totaltime)