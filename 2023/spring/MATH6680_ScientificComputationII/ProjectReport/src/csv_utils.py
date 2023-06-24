import sys
from enum import IntEnum

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

from matplotlib import cm

DATA_FILE = 'data/homogeneous_swe_observation.csv'


class PhysicalQuantity(IntEnum):
    """
    Enum values correspond to array index in loaded CSV data.
    """
    HEIGHT = 0
    VELOCITY = 1
    BATHYMETRY = 2


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
        -sample['alpha'],
    ))  # [(x, t, h, v, alpha)]


def plot_point_cloud(points, quantity=PhysicalQuantity.BATHYMETRY):
    """
    """

    samples = points[1::5]

    if quantity == PhysicalQuantity.BATHYMETRY:
        # Numerically integrate bathymetry tangent angle to get bathymetry curve
        #
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)

        #ax.set_xlabel('x')
        #ax.set_ylabel('B(x)')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$\beta(x)$')
        ax.set_ylim([-0.2, 0.2])

        initial_samples = samples[np.where(np.isclose(samples[:, 1], 0))]
        initial_bathymetry = initial_samples[:, 4]
        xx = initial_samples[:, 0]

        # We plot against the negative bathymetry angle, since our training data came from
        # alpha measured in the opposite direction.
        #
        ax.plot(
            xx,
            -initial_bathymetry,
            # Optional: Numerically integrate to get bathymetry
            #
            #xx[:-1],
            #scipy.integrate.cumulative_trapezoid(np.tan(-initial_bathymetry), x=xx),
            marker='o',
        )
    else:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        ax.set_xlabel('x')
        ax.set_ylabel('t')

        u = samples[:, quantity+2]

        ax.plot_trisurf(
            samples[:, 0],
            samples[:, 1],
            u,
            cmap=cm.viridis,
            linewidth=0,
            antialiased=False
        )
        m = cm.ScalarMappable(cmap=cm.viridis)
        m.set_array(u)

        plt.colorbar(m)
    
    if quantity == PhysicalQuantity.BATHYMETRY:
        plt.savefig('swe_pseudospectral_bathymetry.png', bbox_inches='tight')
    elif quantity == PhysicalQuantity.HEIGHT:
        plt.savefig('swe_pseudospectral_height.png', bbox_inches='tight')
    elif quantity == PhysicalQuantity.VELOCITY:
        plt.savefig('swe_pseudospectral_velocity.png', bbox_inches='tight')

    plt.show()


if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) > 1:
        plot_point_cloud(load_data(file=args[1])[1::5], quantity=int(args[0], 10))
    elif len(args) > 0:
        plot_point_cloud(load_data()[1::5], quantity=int(args[0], 10))
    else:
        plot_point_cloud(load_data()[1::5])
