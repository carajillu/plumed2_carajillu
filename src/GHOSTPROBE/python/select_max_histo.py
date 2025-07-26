import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def histogram_peak_bin_center(df: pd.DataFrame, column: str, nbins: int = 10, plot: bool = True) -> float:
    """
    Generate a histogram of values from a DataFrame column, and return the center
    of the bin with the highest number of elements.

    Parameters:
    - df: pandas DataFrame
    - column: name of the column to generate the histogram for
    - nbins: number of bins in the histogram
    - plot: if True, display the histogram

    Returns:
    - The center of the bin with the highest number of elements
    """
    values = df[column].dropna().values
    counts, bin_edges = np.histogram(values, bins=nbins)

    max_bin_index = np.argmax(counts)
    bin_center = (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1]) / 2

    if plot:
        plt.hist(values, bins=nbins, edgecolor='black')
        plt.axvline(bin_center, color='red', linestyle='--', label=f'Peak bin center: {bin_center:.2f}')
        plt.xlabel(column)
        plt.ylabel('Count')
        plt.title(f'Histogram of {column}')
        plt.legend()
        plt.show()

    return bin_center
