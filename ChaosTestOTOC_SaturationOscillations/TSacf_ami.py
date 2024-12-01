import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mutual_info_score
from sklearn.preprocessing import KBinsDiscretizer
from statsmodels.tsa.stattools import acf
import pandas as pd
import os

# Function to compute average_mutual_information using histograms
def average_mutual_information(time_series, max_lag, n_bins=30):
    n = len(time_series)
    ami_values = np.zeros(max_lag)
    discretizer = KBinsDiscretizer(n_bins=n_bins, encode='ordinal', strategy='uniform')
    time_series_binned = discretizer.fit_transform(time_series.reshape(-1, 1)).flatten()
    
    for lag in range(1, max_lag + 1):
        x = time_series_binned[:-lag]
        y = time_series_binned[lag:]
        ami_values[lag - 1] = mutual_info_score(x, y)
        
    return ami_values

# Load time series data from file

mu = 20


for temp in range(2,11,2):
        
           
    filename = f'newCT200Bmu_{mu}_T={temp}.csv'
    
    
    data = pd.read_csv(filename, delim_whitespace=True, header=None)
    
    end = None
    t = data[0].str.split(',').str[0].astype(float).to_numpy()
    OTOC_values = data[0].str.split(',').str[1].astype(float).to_numpy()
    
    # Compute average mutual information
    max_lag = 2000
    ami_values = average_mutual_information(OTOC_values, max_lag)
    # Compute auto-correlation function
    acf_values = acf(OTOC_values, nlags=max_lag, fft=True)


    plt.suptitle(f'$\mu$ ={mu}, $T$ = {temp}')
    plt.subplot(1,2,1)
    plt.plot(range(1,max_lag+1),ami_values, color='blue')
    plt.xlabel('Lag')
    plt.ylabel('AMI')
    
    #plt.subplot(1,2,1)
    #plt.plot(t[0:390], OTOC_values[0:390])
    
    plt.subplot(1,2,2)
    plt.plot(range(1,max_lag+2),acf_values, color='red')
    plt.tight_layout()
    plt.xlabel('lag')
    plt.ylabel('ACF')
    

    folder_path = 'C:\\Users\\heckm\\OneDrive\\Documents\\Python Scripts\\Project 1\\Time Series\\OTOC_results'
    # Make sure the folder exists
    os.makedirs(folder_path, exist_ok=True)
    
    # Define the file name
    file_name = f'AMI_ACF_mu_{mu}_T={temp}.png'
    # Combine folder path and file name
    file_path = os.path.join(folder_path, file_name)
    # Save the figure
    plt.savefig(file_path)
    plt.show()
