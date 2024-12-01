import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
from scipy.fft import fft, fftfreq
from scipy.signal import welch
import os

mu = 0


for temp in range(2,11,2):
             
    filename = f'newCT200Bmu_{mu}_T={temp}.csv'
    
    #os.mkdir(f'C:\\Users\\heckm\\OneDrive\\Documents\\Python Scripts\\Project 1\\Time Series\\OTOC_results\\LLE_fig_mu_{mu}_T_{temp}')
    data = pd.read_csv(filename, sep='\s+', header=None)
    time = data[0].str.split(',').str[0].astype(float)
    OTOC_values = data[0].str.split(',').str[1].astype(float)
    
    t = time.to_numpy()
    time_series = OTOC_values.to_numpy()
    end = None
    
    
    #Shows a normal takens' embedding with tau as number of data points.
    def embed_time_series(time_series, dim, tau):
        n = len(time_series)
        if n - (dim - 1) * tau <= 0:
            raise ValueError("Time series is too short for given dim and tau")
        return np.array([time_series[i:i + (dim - 1) * tau + 1:tau] for i in range(n - (dim - 1) * tau)])
    
    delay = 3
    embedded=embed_time_series(time_series[1:end],8,delay)
    
    # Compute FFT and power spectrum
    fft_values = fft(time_series)
    
    power_spectrum = np.abs(fft_values)**2 / len(time_series)
    
    freqs = fftfreq(len(time_series), d=(t[1] - t[0]))
    
    positive_freqs = freqs[:len(freqs)//2]
    positive_power_spectrum = power_spectrum[:len(power_spectrum)//2]
    
    mean_frequency = np.sum(positive_freqs * positive_power_spectrum) / np.sum(positive_power_spectrum)
    
    print("Mean Frequency:", mean_frequency) 
    min_dist = np.ceil((1/(mean_frequency * (1/(t[1]-t[0]))))).astype(int)
    if min_dist < 10:
        min_dist = 10
        
        
    max_iter = 1000
    iteration_count =-1   
    timespan = 500
    LLEs = np.zeros((timespan,len(range(0,10000,1000))))
    coefficients = np.zeros((2,len(range(0,10000,1000))))
    
    
    for move in range(0,10000,1000):
        
        iteration_count += 1
            
        distcheck_euc_extended = cdist(embedded[0+move:2500+move], embedded[0+move:2500+move], 'euclidean')    
        np.fill_diagonal(distcheck_euc_extended, np.inf) 
        distMat_mindist = np.zeros((max_iter-min_dist,max_iter))
        
        for i in range(len(distMat_mindist)):
            distMat_mindist[i] = distcheck_euc_extended[i,i+min_dist:(i+min_dist+max_iter)]
        
        neighbors_euc_test = np.argmin(distMat_mindist, axis=1)
        neighbors_euc_vals_test=np.min(distMat_mindist,axis=1)
        
        for i in range(len(neighbors_euc_test)):
            neighbors_euc_test[i] = neighbors_euc_test[i]+i+min_dist
            
       
        
            ####
        def joes_lyapunov(time,max_iter,inv_sampling_freq):
            joes_divsum= np.zeros((time,max_iter))
            for j in range(max_iter):
                print(j)
                for t in range(time):
                    print(t)
                    distance_after_timesteps_t=distcheck_euc_extended[j+t,neighbors_euc_test[j]+t]
                    print(j+t)
                    print(neighbors_euc_test[j]+t)
                    joes_divsum[t,j] = np.log(distance_after_timesteps_t) ##ln of distance between nearest neighbors at (time)/ their initial distance
        
        
            joes_divsum_rows = np.sum(joes_divsum, axis= 1)
            joes_divsum_final = joes_divsum_rows / (max_iter)
            return joes_divsum, joes_divsum_rows, joes_divsum_final
        
        
        test, test_rows_sum, test_final = joes_lyapunov(timespan,max_iter-min_dist,(t[1]-t[0]))
        
        
    
        folder_path = f'C:\\Users\\heckm\\OneDrive\\Documents\\Python Scripts\\Project 1\\Time Series\\OTOC_results\\LLE_fig_mu_{mu}_T_{temp}'
        file_name = f'LLEfig_mu_{mu}_T={temp}_{iteration_count}.png'
        file_path = os.path.join(folder_path, file_name)
        os.makedirs(folder_path, exist_ok=True)
        
        t = t.astype(float)[0:timespan]
        plt.rcParams['font.family'] = 'Times New Roman'
        plt.plot(t,test_final, color='purple')
        plt.xlabel('Time (S)', fontsize=10)
        plt.ylabel('Log(divergence), Ensemble Average',fontsize=10)
        plt.title(f'$\mu$ = {mu}, $T$ = {temp}, delay = {delay}')
        plt.xticks(fontsize=8)
        plt.grid(True)
        plt.yticks(fontsize=8)
        #plt.savefig(file_path)
        
        
        LLEs[:,iteration_count] = test_final
        plt.show()
    
        #coefficients[:,iteration_count] = np.polyfit(t[0:50], test_final[0:50],1)
    LLEs_df = np.savetxt(f'LLEs_mu_{mu}_T_{temp}_numIter_{max_iter}.dat',LLEs,delimiter=',')
