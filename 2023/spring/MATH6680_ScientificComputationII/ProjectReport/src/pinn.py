import deepxde as dde
import numpy as np

from deepxde.backend import tf

def swe_1d(domain, ics, bcs, g=9.81, C_f=0, samples=None, iterations=50000, lr=10e-3, horizontal=False):
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
            (v_t + (v * v_x + g * h_x * tf.sin(alpha + np.pi/2))) - (g * tf.sin(alpha) * (1 + .5 * h * alpha_x) - C_f * v * v / h),
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
    
        # Wave height network architecture [20 x 3]
        #
        h_layer_size = 20 
    
        h = tf.layers.dense(x_t, h_layer_size, tf.nn.tanh)
        h = tf.layers.dense(h, h_layer_size, tf.nn.tanh)
        h = tf.layers.dense(h, h_layer_size, tf.nn.tanh)
        h_output = tf.layers.dense(h, 1, None)
    
        # Velocity network architecture [20 x 3]
        #
        v_layer_size = 20
    
        v = tf.layers.dense(x_t, v_layer_size, tf.nn.tanh)
        v = tf.layers.dense(v, v_layer_size, tf.nn.tanh)
        v = tf.layers.dense(v, v_layer_size, tf.nn.tanh)
        v_output = tf.layers.dense(v, 1, None)
    
        # Bathymetry angle network architecture [20 x 3]
        #
        alpha_layer_size = 20
    
        alpha = tf.layers.dense(x, alpha_layer_size, tf.nn.tanh)
        alpha = tf.layers.dense(alpha, alpha_layer_size, tf.nn.tanh)
        alpha = tf.layers.dense(alpha, alpha_layer_size, tf.nn.tanh)
        alpha_output = tf.layers.dense(alpha, 1, None)
    
        final_output = tf.concat([
            h_output, 
            v_output, 
            x*0 if horizontal else alpha_output,
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

    # Predict neural network solution on random points and display residual value
    #     
    X = domain.random_points(100000)
 
    f = model.predict(X, operator=_pde)
    err_eq = np.absolute(f)
    err = np.mean(err_eq)
    print('Mean residual: %.3e' % (err))
    
    dde.saveplot(losshistory, train_state, issave=True, isplot=True)

    # Compute L2 error
    #  
    if samples is not None and len(samples):
        y_true = samples[:, 2:5]
        y_pred = model.predict(samples[:, 0:2])
        print('L2 relative error:', dde.metrics.l2_relative_error(y_true, y_pred))

    return model

# Example usage: solves the 1D SWE with horizontal bathymetry, periodic boundary conditions, and initial conditions:
#
# u(x, 0) = 2 + sin(pi*x/100)
# v(x, 0) = 0.
#
if __name__ == '__main__':
    # Set random seed to 0 to allow for reproducable randomness in results
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
    bc_h = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=0, derivative_order=0)
    bc_v = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=1, derivative_order=0)

    # Periodic boundary conditions for height and velocity derivatives, respectively
    #
    bc_h_x = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=0, derivative_order=1)
    bc_v_x = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=1, derivative_order=1)

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
        iterations=20000,
        horizontal=True,
    )
