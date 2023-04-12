from enum import Enum

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

DATA_FILE = 'data/linearBsineSWE.csv'

class PhysicalQuantity(Enum):
    HEIGHT = 1
    VELOCITY = 2
    BATHYMETRY = 3

def load_data(sample_size=None, file=DATA_FILE):
    data = np.genfromtxt(file, dtype=float, delimiter=',', names=True)

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

def plot_point_cloud(points, quantity=PhysicalQuantity.HEIGHT):
    """
    """
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    samples = points[1::5]

    if quantity == PhysicalQuantity.HEIGHT:
        u = samples[:, 2:3]
    elif quantity == PhysicalQuantity.VELOCITY:
        u = samples[:, 3:4]
    elif quantity == PhysicalQuantity.BATHYMETRY:
        u = samples[:, 4:5]
    else:
        u = samples[:, 2:3]


    ax.scatter(
        samples[:, 0:1], 
        samples[:, 1:2], 
        u,
        marker='o',
    )

    plt.show()


if __name__ == '__main__':
    plot_point_cloud(load_data())
