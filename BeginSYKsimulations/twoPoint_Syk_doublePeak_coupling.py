##first deformation from normal coupling distribution to a bimodal one. 
##Motivation for this is that some n->n+1 ads/cft might recreate SYK but with 'hypergeometric' coupling distirbution. 
##However, a continuous analog to hypergeom isn't totally acessible. 
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import time
import random
start_time = time.time()
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


Hamiltonian = np.zeros((int(op_size),int(op_size)),dtype=complex)
std_dev = np.sqrt(((6))/(N_maj**3))

mu1 = -2.5*std_dev
mu2 = 2.5*std_dev

##make arrays of fermions 
Fermions = np.vstack((chiEven_loop,chiOdd_loop))


##classic 4 interaction case. can use a gigantic computer to simulate a huge number of interaction terms with very large operators. 
for i in range(0,int(N_maj)):
    for j in range(0,int(N_maj)):
        for k in range(0,int(N_maj)):
            for l in range(0,int(N_maj)):
                print(f'Chi{i}*Chi{j}*Chi{k}*Chi{l}')
                #randomNumber = np.random.normal(mean,std_dev)
                rand = random.randint(1,2)
                
                if rand ==1:
                    Hamiltonian[:,:] += np.random.normal(mu1,std_dev)*(Fermions[i,:,:] @ Fermions[j,:,:] @ Fermions[k,:,:] @ Fermions[l,:,:])
                if rand ==2:
                    Hamiltonian[:,:] += np.random.normal(mu2,std_dev)*(Fermions[i,:,:] @ Fermions[j,:,:] @ Fermions[k,:,:] @ Fermions[l,:,:])
                    
                    
                #Hamiltonian[:,:] += randomNumber*(Fermions[i,:,:] @ Fermions[j,:,:] @ Fermions[k,:,:] @ Fermions[l,:,:])
                
##this is for 4 interaction term, (i^(interaction terms/2))
Hamiltonian = (-1/24)*Hamiltonian
# diving into the 2point here 
beta = 1
range_values = np.arange(0,beta,.01) ##Periodicity at beta right?
array_length = len(range_values)
new_array2 = np.zeros((int(N_maj),array_length))
new_array2Test = np.zeros((int(N_maj),array_length))

imaginary_time_evo = linalg.expm(-beta*Hamiltonian) 
for i in range(0,int(N_maj)):   
    count = -1
    for tau in np.arange(0,beta,.01):
        count += 1
        
        #time_evo_operator is its own thing to keep all of the time evolution parts attached to a Majorana operator
        time_evo_operator = linalg.expm((-beta+tau)*Hamiltonian) @ Fermions[i,:,:] @ linalg.expm(-1*Hamiltonian*tau)  ##they win again, that being the youtube video Leo and Evan keep telling me to watch and GPT
        point = np.trace(time_evo_operator @ Fermions[i,:,:]) / np.trace(linalg.expm((-beta+tau)*Hamiltonian) @ linalg.expm(-1*Hamiltonian*tau))
            
    
        new_array2[i,count]=point

twoPoint = np.mean(new_array2,axis=0)

##plot what you got 
plt.rc('font', family='Times New Roman')
plt.plot(range_values,twoPoint, color='purple')
plt.title(r'$G_{R} $', fontsize = 14)
plt.xlabel(r'$\tau$', fontsize= 13)
plt.ylabel(r'Two-Point($\tau$)', fontsize = 13)
plt.grid()
plt.show()

end_time = time.time()
total_time = end_time-start_time
print(f"Time for {N_maj} Fermions: {total_time} seconds")

##now save 
##shows perioditcity through tau = 0 to tau = beta. Means something. 