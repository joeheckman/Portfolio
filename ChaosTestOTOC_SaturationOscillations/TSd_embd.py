import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

def embed_time_series(time_series, dim, tau):
    n = len(time_series)
    if n - (dim - 1) * tau <= 0:
        raise ValueError("Time series is too short for given dim and tau")
    return np.array([time_series[i:i + (dim - 1) * tau + 1:tau] for i in range(n - (dim - 1) * tau)])

def false_nearest_neighbors(time_series, dim, tau, rtol=15.0, atol=2.0):
    n = len(time_series)
    fnn_percentages = []
    
    embedded = embed_time_series(time_series, dim, tau)
    embedded_next = embed_time_series(time_series, dim + 1, tau)

    # Ensure both arrays have compatible shapes
    if embedded_next.shape[0] < n - (dim + 1 - 1) * tau:
        return fnn_percentages  # Return empty if not enough points

    distances = np.linalg.norm(embedded[:, None] - embedded, axis=2)
    np.fill_diagonal(distances, np.inf)  # avoid self-neighbors
    nearest_indices = np.argmin(distances, axis=1)

    # Ensure nearest_indices are within bounds of embedded_next
    valid_indices = np.arange(embedded_next.shape[0])
    nearest_indices = np.clip(nearest_indices, 0, valid_indices[-1])

    nearest_distances = distances[np.arange(n - (dim - 1) * tau), nearest_indices]
    next_distances = np.linalg.norm(embedded_next[-1] - embedded_next[nearest_indices], axis=1)

    # Handle division by zero or small values
    epsilon = 1e-10  # Small value to avoid division by zero
    valid_mask = nearest_distances > epsilon
    if np.any(valid_mask):
        fnn = (next_distances[valid_mask] / nearest_distances[valid_mask] > rtol) | (next_distances[valid_mask] > atol)
        fnn_percent = np.sum(fnn) / len(fnn) * 100
    else:
        fnn_percent = np.nan

    fnn_percentages.append(fnn_percent)

    return fnn_percentages

# Parameters
max_dim = 10
max_tau = 10
rtol = 15.0
atol = 2.0
max_length = 5000

# Load time series data from file
filename = 'TSlorenz.dat'  # Replace with your file name
data = pd.read_csv(filename, sep='\s+', header=None)
t = data.iloc[1:4999, 0].values
time_series = np.array(data.iloc[1:4999, 1].values,dtype = float)

# Check length of time series and interpolate if necessary
if len(time_series) > max_length:
    interpolation_function = interp1d(t, time_series, kind='cubic')
    new_t = np.linspace(t.min(), t.max(), max_length)
    time_series = interpolation_function(new_t)
    t = new_t

# Prepare data for 2D plot
dims = np.arange(1, max_dim + 1)
taus = np.arange(1, max_tau + 1)
fnn_results = []

for dim in dims:
    fnn_for_dim = []
    for tau in taus:
        fnn_percentages = false_nearest_neighbors(time_series, dim, tau, rtol, atol)
        if fnn_percentages:
            fnn_for_dim.append(fnn_percentages[0])
        else:
            fnn_for_dim.append(np.nan)  # Use NaN for missing values
    fnn_results.append(fnn_for_dim)

fnn_results = np.array(fnn_results)

# Create 2D plot
plt.figure(figsize=(10, 6))
for dim in range(max_dim):
    plt.plot(taus, fnn_results[dim], label=f'dim={dim+1}')

plt.xlabel('Time Delay (tau)')
plt.ylabel('FNN Percentage')
plt.title('False Nearest Neighbors Analysis')
plt.legend()
plt.grid(True)
plt.show()

