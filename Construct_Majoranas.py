"""
Create Majoranas
"""

#guessing here boys easy does it
import numpy as np

N_maj = 20 ##trivially simple (haha dick)

##find dtype=complex is the homework problem
op_size = 2**(N_maj/2)
print([f'size of operator is {op_size} x {op_size}'])
chiEven = np.zeros((int(N_maj/2),int(op_size),int(op_size)),dtype=complex)
chiOdd = np.zeros((int(N_maj/2),int(op_size),int(op_size)),dtype=complex)
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