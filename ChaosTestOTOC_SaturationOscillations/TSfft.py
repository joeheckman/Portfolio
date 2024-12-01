import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import welch
import pandas as pd
import os 




# Load time series data from file
#filename = 'TSunk.dat'  # Replace with your file name
#data = pd.read_csv(filename, sep='\s+', header=None)
#end = None
#t = np.array(data.iloc[1:end, 0].values, dtype=float)
#signal = data.iloc[1:end, 1].values


for temp in range(2,11,2):
    ##for the .csv files 
    filename = f'newCT200Bmu_20_T={temp}.csv'
    data = pd.read_csv(filename, sep='\s+', header=None)
    time = data[0].str.split(',').str[0].astype(float)
    OTOC_values = data[0].str.split(',').str[1].astype(float)
    
    t = time.to_numpy()
    signal = OTOC_values.to_numpy()
    end = None
    
    
    # Compute FFT and power spectrum
    fft_values = fft(signal)
    
    power_spectrum = np.abs(fft_values)**2 / len(signal)
    
    freqs = fftfreq(len(signal), d=(t[1] - t[0]))
    
    positive_freqs = freqs[:len(freqs)//2]
    positive_power_spectrum = power_spectrum[:len(power_spectrum)//2]
        
        
    # Plot the power spectrum
    plt.figure(figsize=(10, 6))
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.plot(positive_freqs[1:end], positive_power_spectrum[1:end], color='purple')
    plt.title(f'Power Spectrum of Time Series $\mu$ = 20, $T$ = {temp}')
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    #plt.xlim(0, 0.5 / (t[1] - t[0]))  # Limit the x-axis to focus on lower frequencies
    plt.xlim(0, .1)
    plt.grid(True)
    
    
    folder_path = 'C:\\Users\\heckm\\OneDrive\\Documents\\Python Scripts\\Project 1\\Time Series\\OTOC_results\\power_spectrums'
    
    # Make sure the folder exists
    os.makedirs(folder_path, exist_ok=True)
    
    # Define the file name
    file_name = f'newCT200Bmu_20_T={temp}_fft.png'
    # Combine folder path and file name
    file_path = os.path.join(folder_path, file_name)
    # Save the figure
    plt.savefig(file_path)
    plt.show()
