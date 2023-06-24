import math
import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from matplotlib import cm
from deepxde.backend import tf

from csv_utils import load_data, PhysicalQuantity

print(tf.config.list_physical_devices('GPU'))


def swe_1d(domain, ics, bcs, g=9.81, C_f=0, samples=None, iterations=50000, lr=10e-3,
           prescribed_bathymetry=None):
    def _pde(x, y):
        """
        Encodes the SWE residual in terms of output neurons (i.e., the neurons in the final network layer which
        correspond to h, v, and alpha).
        """
        h = y[:, 0:1]      # wave height
        v = y[:, 1:2]      # velocity
        alpha = y[:, 2:3]  # bathymetry

        h_x = dde.grad.jacobian(y, x, i=0, j=0)
        h_t = dde.grad.jacobian(y, x, i=0, j=1)

        v_x = dde.grad.jacobian(y, x, i=1, j=0)
        v_t = dde.grad.jacobian(y, x, i=1, j=1)

        alpha_x = dde.grad.jacobian(y, x, i=2, j=0)

        # Inhomogeneous system with variable bathymetry (alpha)
        #
        return [
            (h_t + (h * v_x + h_x * v)) + 0,
            (v_t + (v * v_x + g * h_x * tf.math.cos(alpha) - g * h * tf.math.sin(alpha) * alpha_x))
            - g * tf.math.sin(alpha)
            + C_f * v * v / h
            + 0.5 * g * h * tf.math.sin(alpha) * alpha_x,
        ]

    def _modify_network_output(input, *args, **kwargs):
        """
        We implement a custom neural network architecture in this function so that we can restrict the bathymetry
        function to a spatial domain since bathymetry (presumably) does not vary with time.
        """
        x = input[:, 0:1]
        time = input[:, 1:2]

        x_t = tf.concat([
            x,
            time,
        ], axis=1)

        # Wave height/velocity network architecture [20 x 4]
        #
        hv_layer_size = 20

        hv = tf.layers.dense(x_t, hv_layer_size, tf.nn.tanh)
        hv = tf.layers.dense(hv, hv_layer_size, tf.nn.tanh)
        hv = tf.layers.dense(hv, hv_layer_size, tf.nn.tanh)
        hv_output = tf.layers.dense(hv, 2, None)

        # Bathymetry angle network architecture [20 x 4]
        #
        alpha_layer_size = 20

        alpha = tf.layers.dense(x, alpha_layer_size, tf.nn.tanh)
        alpha = tf.layers.dense(alpha, alpha_layer_size, tf.nn.tanh)
        alpha = tf.layers.dense(alpha, alpha_layer_size, tf.nn.tanh)
        alpha_output = tf.layers.dense(alpha, 1, None)

        final_output = tf.concat([
            hv_output[:, 0:1],
            hv_output[:, 1:2],
            alpha_output if prescribed_bathymetry is None else prescribed_bathymetry(x),
        ], axis=1)

        return final_output

    # Encodes ICs and BCs into a loss function to train the neural network solution
    #
    data = dde.data.TimePDE(
        domain,
        _pde,
        [
          *ics,
          *bcs,
        ],
        num_domain=5000,
        num_boundary=200,
        num_initial=320,
    )

    # Neural network structure
    #
    network_input_outputs = [2, 3]

    # Create the neural network solution surrogate
    #
    net = dde.nn.FNN(
        network_input_outputs,
        'tanh',                # Activation function
        'Glorot normal',       # Weight initialization scheme
    )

    net.apply_output_transform(_modify_network_output)
    model = dde.Model(data, net)

    # Train network with ADAM and L-BFGS optimizers
    #
    model.compile('adam', lr=lr)
    model.train(iterations=iterations)
    model.compile('L-BFGS')
    losshistory, train_state = model.train()

    # Predict neural network solution on random points (if samples are not provided) and display
    # residual value
    #
    X = samples[:, 0:2] if samples is not None else domain.random_points(100000)

    f = model.predict(X, operator=_pde)
    err_eq = np.absolute(f)
    err = np.mean(err_eq)
    print('Mean residual: %.3e' % (err))

    dde.saveplot(losshistory, train_state, issave=True, isplot=True)

    # Compute L2 error if sample data is provided
    #
    if samples is not None and len(samples):
        y_true = samples[:, 2:5]
        y_pred = model.predict(samples[:, 0:2])

        print('L2 relative error:', dde.metrics.l2_relative_error(y_true, y_pred))
        print('Plots at sampled points:')

        for i in (0, 1, 2):
            if i == PhysicalQuantity.BATHYMETRY:
                # Numerically integrate bathymetry tangent angle to get bathymetry curve
                #
                fig = plt.figure(figsize=(8, 8))
                ax = fig.add_subplot(111)

                ax.set_xlabel(r'$x$')
                ax.set_ylabel(r'$\beta(x)$')
                ax.set_ylim([-0.2, 0.2])

                initial_bathymetry_idx = np.isclose(samples[:, 1], 0)
                bathymetry_angle = y_pred[np.where(initial_bathymetry_idx)][:, PhysicalQuantity.BATHYMETRY]
                xx = samples[np.where(initial_bathymetry_idx)][:, 0]

                # We plot against the negative bathymetry angle, since our training data came from
                # alpha measured in the opposite direction.
                #
                ax.plot(
                    xx,
                    -bathymetry_angle,
                    # Optional: Numerically integrate to get bathymetry
                    #
                    #xx[:-1],
                    #scipy.integrate.cumulative_trapezoid(np.tan(-bathymetry_angle), x=xx),

                    marker='o',
                )
            else:
                fig = plt.figure(figsize=(8, 8))
                ax = fig.add_subplot(111, projection='3d')

                ax.set_xlabel('x')
                ax.set_ylabel('t')
                # ax.set_zlim([-np.pi, np.pi])

                ax.plot_trisurf(
                    samples[:, 0],
                    samples[:, 1],
                    y_pred[:, i],
                    cmap=cm.viridis,
                    linewidth=0,
                    antialiased=False
                )
                m = cm.ScalarMappable(cmap=cm.viridis)
                m.set_array(y_pred[:, i])
                fig.colorbar(m)

            if i == PhysicalQuantity.BATHYMETRY:
                plt.savefig('swe_pinn_bathymetry.png', bbox_inches='tight')
            elif i == PhysicalQuantity.HEIGHT:
                plt.savefig('swe_pinn_height.png', bbox_inches='tight')
            elif i == PhysicalQuantity.VELOCITY:
                plt.savefig('swe_pinn_velocity.png', bbox_inches='tight')
            plt.show()

    return model


# Example usage: solves the 1D SWE with horizontal bathymetry, periodic boundary conditions, and
# initial conditions:
#
# u(x, 0) = 2 + sin(pi*x/100)
# v(x, 0) = 0.
#
if __name__ == '__main__':
    # Set random seed to 0 to allow for reproducible randomness in results
    #
    dde.config.set_random_seed(0)

    # Set default float to 32-bit since Tensorflow seems to like that better
    # when training on a GPU
    #
    dde.config.set_default_float('float32')

    # Generate domain
    #
    x_min, x_max = (0, 200)
    t_min, t_max = (0, 15)

    geom = dde.geometry.Interval(x_min, x_max)
    timedomain = dde.geometry.TimeDomain(t_min, t_max)
    geomtime = dde.geometry.GeometryXTime(geom, timedomain)

    # Periodic boundary conditions for height and velocity, respectively
    #
    bc_h = dde.icbc.PeriodicBC(
        geomtime,
        0,
        lambda _, on_boundary: on_boundary,
        component=0,
        derivative_order=0
    )
    bc_v = dde.icbc.PeriodicBC(
        geomtime,
        0,
        lambda _, on_boundary: on_boundary,
        component=1,
        derivative_order=0,
    )

    # Periodic boundary conditions for height and velocity derivatives, respectively
    #
    bc_h_x = dde.icbc.PeriodicBC(
        geomtime,
        0,
        lambda _, on_boundary: on_boundary,
        component=0,
        derivative_order=1,
    )
    bc_v_x = dde.icbc.PeriodicBC(
        geomtime,
        0,
        lambda _, on_boundary: on_boundary,
        component=1,
        derivative_order=1,
    )

    # Initial condition for wave height
    #
    ic_h = dde.icbc.IC(
        geomtime,
        lambda x: 2 + np.sin(x[:, 0:1] * np.pi / 100),
        lambda _, on_initial: on_initial,
        component=0,
    )

    # Initial condition for wave velocity
    #
    ic_v = dde.icbc.IC(
        geomtime,
        lambda x: 0,
        lambda _, on_initial: on_initial,
        component=1,
    )

    # Solve system
    #
    swe_1d(
        geomtime,
        [
            ic_h,
            ic_v
        ],
        [
            bc_h,
            bc_v,
            bc_h_x,
            bc_v_x,
        ],
        iterations=10000,
        prescribed_bathymetry=lambda x: x*0,  # horizontal bathymetry
    )
