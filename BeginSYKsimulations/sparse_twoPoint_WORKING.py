import numpy as np
import scipy.sparse
import scipy.linalg as linalg
import matplotlib.pyplot as plt



N_maj = 10


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
        temp = np.kron(temp,chi_even_builder[i+1+count]); ##this satisfies N_maj = 6
        temp_odd = np.kron(temp_odd,chi_odd_builder[i+1+count])
        print(i+1+count)
        
        if len(temp) == op_size:
            print(op_size)
            chiEven_loop[mat_pos,:,:] = temp;
            chiOdd_loop[mat_pos,:,:] = temp_odd;
            break


##test clifford algebras
clif_test_1 = np.zeros((int(N_maj/2),int(op_size),int(op_size)),dtype=complex)
clif_test_2 = np.zeros((int(N_maj/2),int(N_maj/2),int(op_size),int(op_size)),dtype=complex)
clif_test_1_odd = np.zeros((int(N_maj/2),int(op_size),int(op_size)),dtype=complex)
clif_test_2_odd = np.zeros((int(N_maj/2),int(N_maj/2),int(op_size),int(op_size)),dtype=complex)

for test in range(0,int(N_maj/2)):
    clif_test_1[test,:,:] = np.dot(chiEven_loop[test,:,:],chiEven_loop[test,:,:])+np.dot(chiEven_loop[test,:,:],chiEven_loop[test,:,:]) ##test for correct construction
    clif_test_1_odd[test,:,:] = np.dot(chiOdd_loop[test,:,:],chiOdd_loop[test,:,:])+np.dot(chiOdd_loop[test,:,:],chiOdd_loop[test,:,:])
##hold on
for p in range(0,int(N_maj/2)):
    for q in range(0,int(N_maj/2)):
        clif_test_2[p,q,:,:] = np.dot(chiEven_loop[p,:,:],chiEven_loop[q,:,:])+np.dot(chiEven_loop[q,:,:],chiEven_loop[p,:,:]) ##test for correct construction
        ##clif_test_2_odd[p,q,:,:] satisfies clif_test_1 when p = q, therefore the 2$\delta_{ij}$
        clif_test_2_odd[p,q,:,:] = np.dot(chiOdd_loop[p,:,:],chiOdd_loop[q,:,:])+np.dot(chiOdd_loop[q,:,:],chiOdd_loop[p,:,:])
        
##now let's make the operators into 'Fermions' because that makes more sense in my head and I think that's what other people do too.

Hamiltonian = np.zeros((int(op_size),int(op_size)),dtype=complex)
HamiltonianTest = np.zeros((int(op_size),int(op_size)),dtype=complex)
mean = 0
std_dev = np.sqrt(((6))/(16*N_maj**3))

##make sparse arrays of fermions 
Fermions = np.vstack((chiEven_loop,chiOdd_loop))
sparseFermions = [scipy.sparse.csr_matrix(Fermions[i]) for  i in range(N_maj)]


##classic 4 interaction case. can use a gigantic computer to simulate a huge number of interaction terms with very large operators. 
for i in range(0,int(N_maj)):
    for j in range(0,int(N_maj-1)):
        for k in range(0,int(N_maj-2)):
            for l in range(0,int(N_maj-3)):
                print(f'Chi{i}*Chi{j}*Chi{k}*Chi{l}')
                randomNumber = np.random.normal(mean,std_dev)
                
                HamiltonianTest[:,:] += randomNumber*(Fermions[i,:,:] @ Fermions[j,:,:] @ Fermions[k,:,:] @ Fermions[l,:,:])
                Hamiltonian[:,:] += randomNumber*(sparseFermions[i] @ sparseFermions[j] @ sparseFermions[k] @ sparseFermions[l])

##this is for 4 interaction term, (i^(interaction terms/2))
Hamiltonian = -1*Hamiltonian
#2point, unnormalized but it's only a factor 
beta = 1
range_values = np.arange(0,beta,.01) ##Periodicity at beta right?
array_length = len(range_values)
new_array2 = np.zeros((int(N_maj),array_length,int(N_maj)))
new_array2Test = np.zeros((int(N_maj),array_length,int(N_maj)))

imaginary_time_evo = linalg.expm(-beta*Hamiltonian) 
for i in range(0,int(N_maj)):
    for j in range(0,int(N_maj)):   
        count = -1
        for tau in np.arange(0,beta,.01):
            count += 1
            Testtime_evo_operator = linalg.expm((-beta+tau)*HamiltonianTest) @ Fermions[j,:,:] @ linalg.expm(-1*HamiltonianTest*tau)  ##they win again, that being the youtube video Leo and Evan keep telling me to watch and GPT
            Testpoint = np.trace(Testtime_evo_operator @ Fermions[i,:,:]) / np.trace(linalg.expm((-beta+tau)*HamiltonianTest) @ linalg.expm(-1*HamiltonianTest*tau))
            
            time_evo_operator = linalg.expm((-beta+tau)*Hamiltonian) @ sparseFermions[j] @ linalg.expm(-1*Hamiltonian*tau)  ##they win again, that being the youtube video Leo and Evan keep telling me to watch and GPT
            point = np.trace(time_evo_operator @ sparseFermions[i]) / np.trace(linalg.expm((-beta+tau)*Hamiltonian) @ linalg.expm(-1*Hamiltonian*tau))
            
            
            new_array2[i,count,j]=point
            new_array2Test[i,count,j] = point


##ensemble average 
twoPoint = np.zeros((int(N_maj),array_length))
for mean in range(0,int(N_maj)):        
    twoPoint[mean,:] = np.mean(new_array2[mean,:,:], axis=1)
    
twoPointFinal= np.mean(twoPoint,axis=0)

##plot what you got 
plt.rc('font', family='Times New Roman')
plt.plot(range_values,twoPointFinal, color='purple')
plt.title(r'$G_{R} $', fontsize = 14)
plt.xlabel(r'$\tau$', fontsize= 13)
plt.ylabel(r'Two-Point($\tau$)', fontsize = 13)
plt.grid()
plt.show()

##shows perioditcity through tau = 0 to tau = beta. Means something. 