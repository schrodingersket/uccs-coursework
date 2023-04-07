import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

from deepxde.backend import tf
from csv_utils import load_data

dde.config.set_random_seed(0)
dde.config.set_default_float('float32')

# Physical constants
#
g = 9.81  # Gravitational acceleration
C_f = 0   # Frictional coefficient


x_min, x_max = (0, 200)
t_min, t_max = (0, 15)



def pde(x, y):
    h = y[:, 0:1]      # wave height
    v = y[:, 1:2]      # velocity
    alpha = y[:, 2:3]  # bathymetry

    h_x = dde.grad.jacobian(y, x, i=0, j=0)
    h_t = dde.grad.jacobian(y, x, i=0, j=1)

    v_x = dde.grad.jacobian(y, x, i=1, j=0)
    v_t = dde.grad.jacobian(y, x, i=1, j=1)

    alpha_x = dde.grad.jacobian(y, x, i=2, j=0)

    # Inhomogenous system with variable bathymetry (alpha)
    #
    return [
        (h_t + (h * v_x + h_x * v)) + 0,
        (v_t + (v * v_x + g * h_x * tf.sin(alpha + np.pi/2))) - (g * h * tf.sin(alpha) * (1 - alpha_x) - C_f * v * v / h),
    ]

if __name__ == '__main__':
    geom = dde.geometry.Interval(x_min, x_max)
    timedomain = dde.geometry.TimeDomain(t_min, t_max)
    geomtime = dde.geometry.GeometryXTime(geom, timedomain)
    
    bc_h = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=0, derivative_order=0)
    bc_v = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=1, derivative_order=0)

    bc_h_x = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=0, derivative_order=1)
    bc_v_x = dde.icbc.PeriodicBC(geomtime, 0, lambda _, on_boundary: on_boundary, component=1, derivative_order=1)

    ic_h = dde.icbc.IC(
        geomtime,
        lambda x: 2 + np.sin(x[:, 0:1] * np.pi / 100),
        lambda _, on_initial: on_initial,
        component=0,
    )
    ic_v = dde.icbc.IC(
        geomtime,
        lambda x: 0,
        lambda _, on_initial: on_initial,
        component=1,
    )
    
    data = dde.data.TimePDE(
        geomtime, 
        pde, 
        [
          ic_h, 
          ic_v,
          bc_h,
          bc_v, 
          bc_h_x,
          bc_v_x
        ], 
        num_domain=5000, 
        num_boundary=200, 
        num_initial=320,
    )

    # Neural network structure
    #
    network_input_outputs = [2, 3] 

    # Create neural network solution surrogate
    #
    net = dde.nn.FNN(
        network_input_outputs,
        'tanh',                # Activation function
        'Glorot normal',       # Weight initialization scheme
    )

    def modify_network_output(input, *args, **kwargs):
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
    
        # wave height network architecture [20 x 3]
        #
        h_layer_size = 20 
    
        h = tf.layers.dense(x_t, h_layer_size, tf.nn.tanh)
        h = tf.layers.dense(h, h_layer_size, tf.nn.tanh)
        h = tf.layers.dense(h, h_layer_size, tf.nn.tanh)
        h_output = tf.layers.dense(h, 1, None)
    
        # velocity network architecture [20 x 3]
        #
        v_layer_size = 20
    
        v = tf.layers.dense(x_t, v_layer_size, tf.nn.tanh)
        v = tf.layers.dense(v, v_layer_size, tf.nn.tanh)
        v = tf.layers.dense(v, v_layer_size, tf.nn.tanh)
        v_output = tf.layers.dense(v, 1, None)
    
        # bathymetry angle network architecture [20 x 3]
        #
        alpha_layer_size = 20
    
        alpha = tf.layers.dense(x, alpha_layer_size, tf.nn.tanh)
        alpha = tf.layers.dense(x, alpha_layer_size, tf.nn.tanh)
        alpha = tf.layers.dense(x, alpha_layer_size, tf.nn.tanh)
        alpha_output = tf.layers.dense(alpha, 1, None)
    
        final_output = tf.concat([
            h_output, 
            v_output, 
            x*0, # change to alpha_output when using variable bathymetry
        ], axis=1)
    
        return final_output

    net.apply_output_transform(modify_network_output)
    model = dde.Model(data, net)

    # Train network
    # 
    model.compile('adam', lr=1.0e-3)
    model.train(iterations=50000)
    model.compile('L-BFGS')
    losshistory, train_state = model.train()

    # Predict neural network solution on random points and display residual value
    #     
    X = geomtime.random_points(100000)
 
    f = model.predict(X, operator=pde)
    err_eq = np.absolute(f)
    err = np.mean(err_eq)
    print("Mean residual: %.3e" % (err))
    
    dde.saveplot(losshistory, train_state, issave=True, isplot=True)

    # Compute L2 error
    #  
    samples = load_data(1000)
    y_true = samples[:, 2:4]
    y_pred = model.predict(samples[:, 0:2])
    print('L2 relative error:', dde.metrics.l2_relative_error(samples[:, 2:5], y_pred))
    exit()

    # Plots
    #  
    pred_tmin = 0
    pred_tmax = 1.1
    for time_near in np.arange(pred_tmin, pred_tmax, 0.25):
        predicted_time = np.array(time_near)
    
        for t in sorted(X[:, 1]):
            if time_near >= t:
                predicted_time = t
    
        print('Computing solution at t={}'.format(predicted_time))
        x_pred, t_pred = X[:, 0:1], X[:, 1:2]
        u_pred = y_pred[:, 0:1][t_pred == predicted_time.item()]
        u_obs = y_true[t_pred == predicted_time.item()]
        x_obs = X[:, 0:1][t_pred == predicted_time.item()]
    
        plt.plot(x_obs, u_pred, color='black', label='u-PINN', linewidth=1)
        plt.scatter(x_obs[::4], u_obs[::4], color='red', label='u-obs', marker='x')
        plt.title('u-PINN Predicted IC vs u-Reference IC [t={:0.2f}]'.format(predicted_time.item()), fontsize=9.5)
        plt.legend(loc='lower right')
        plt.xlim(x_min, x_max)
        plt.ylim(-2, 2)
        plt.show()
    