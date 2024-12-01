import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os

def embed_time_series(time_series, dim, tau):
    n = len(time_series)
    if n - (dim - 1) * tau <= 0:
        raise ValueError("Time series is too short for given dim and tau")
    return np.array([time_series[i:i + (dim - 1) * tau + 1:tau] for i in range(n - (dim - 1) * tau)])



def plot_phase_space(time_series, tau, dim):
    embedded = embed_time_series(time_series, dim, tau)
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot only if we have at least 3 dimensions in the embedded data
    if embedded.shape[1] >= 3:
        ax.plot(embedded[0:1000, 0], embedded[0:1000, 1], embedded[0:1000, 2], color='purple', linewidth=.4,)
        ax.set_xlabel('OTOC(t)')
        ax.set_ylabel(f'OTOC(t+{tau})')
        ax.set_zlabel(f'OTOC(t+{2*tau})')
        ax.set_title(f'OTOC Phase Space, Delay={tau}, $\mu$ = {mu}, $T$={temp}')
        
        
        file_name = f'OTOC_mu=20_T={temp}_tau={tau}.png'
        folder_path = 'C:\\Users\\heckm\\OneDrive\\Documents\\Python Scripts\\Project 1\\Time Series\\OTOC_results\\Phase_spaces'
        os.makedirs(folder_path, exist_ok=True)

        # Combine folder path and file name
        file_path = os.path.join(folder_path, file_name)
        # Save the figure
        plt.savefig(file_path)
        
    else:
        raise ValueError("Embedding dimension is too small to plot 3D phase space")


#for mu in range (0,6,5):
mu = 20
    
for temp in range(2,11,2): 
    # Load time series data from file
    filename = f'newCT200Bmu_{mu}_T={temp}.csv'  # Replace with your file name
    #data = pd.read_csv(filename, sep='\s+', header=None, skiprows=1)
    #end = None
    #t = data.iloc[0:3000, 0].values
    #time_series = data.iloc[0:3000, 1].values
    
    
    data = pd.read_csv(filename, sep='\s+', header=None)
    time = data[0].str.split(',').str[0].astype(float)
    OTOC_values = data[0].str.split(',').str[1].astype(float)
    
    t = time.to_numpy()
    time_series = OTOC_values.to_numpy()
    end = None
    
    
    # Parameters
    taus = np.arange(3,12,1)
    embedding_dim = 3
    
    for tau in range(len(taus)):
        plot_phase_space(time_series, taus[tau], embedding_dim)
            