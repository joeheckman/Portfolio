##visualize. See the last figure in 'SYK_Explorations.pdf' for the result.  
import numpy as np
import scipy.linalg as linalg
from scipy.sparse import csr_matrix,vstack,kron, save_npz, load_npz
from scipy.sparse.linalg import expm as sparse_expm
import matplotlib.pyplot as plt

beta = 1
range_values = np.arange(0,beta/2,beta/40) ##Periodicity at beta right?
twoPoint_normal_test = load_npz('1.16twoPointN20.npz')
twoPoint_bimodal_test = load_npz('1.17twoPointN20Bimodal.npz')

twoPoint_normal = twoPoint_normal_test.toarray()
twoPoint_bimodal = twoPoint_bimodal_test.toarray()

twoPointNormal = np.mean(twoPoint_normal,axis=0)
twoPointBimodal = np.mean(twoPoint_bimodal,axis=0)

plt.plot(range_values,twoPointNormal, label='Normal Distribution' , color='purple')
plt.plot(range_values,twoPointBimodal,label='Bimodal Distribution' , color='red')
plt.rc('font', family='Times New Roman')
plt.legend()
plt.legend(loc = 'lower left')
plt.xlabel(r'$\tau$', fontsize= 13)
plt.ylabel(r'Two-Point($\tau$)', fontsize = 13)
plt.grid()
plt.show()
