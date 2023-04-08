import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

DATA_FILE = 'data/homSWESin.csv'

def load_data(sample_size=None):
    data = np.genfromtxt(DATA_FILE, dtype=float, delimiter=',', names=True)

    if sample_size:
        sample = np.random.choice(data, sample_size)
    else:
        sample = data

    return np.column_stack((
        sample['x'], 
        sample['t'], 
        sample['h'], 
        sample['v'], 
        sample['alpha'], 
    ))  # [(x, t, h, v, alpha)]

def plot_point_cloud(points):
    """
    """
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    samples = points[1::5]

    ax.scatter(
        samples[:, 0:1], 
        samples[:, 1:2], 
        samples[:, 4:5],
        marker='o',
    )

    plt.show()


if __name__ == '__main__':
    plot_point_cloud(load_data())
