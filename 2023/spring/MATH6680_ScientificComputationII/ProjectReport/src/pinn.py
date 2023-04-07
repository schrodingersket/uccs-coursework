import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

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
    # alpha = y[:, 2:3]  # bathymetry

    h_x = dde.grad.jacobian(y, x, i=0, j=0)
    h_t = dde.grad.jacobian(y, x, i=0, j=1)

    v_x = dde.grad.jacobian(y, x, i=1, j=0)
    v_t = dde.grad.jacobian(y, x, i=1, j=1)

    # alpha_x = dde.grad.jacobian(y, x, i=2, j=0)

    # Homogenous system with variable bathymetry (alpha)
    #
    return [
        h_t + (h * v_x + h_x * v),
        v_t + (v * v_x + g * h_x), 
    ]

    # Inhomogenous system with variable bathymetry (alpha)
    #
    return [
        (h_t + (h * v_x + h_x * v)) + 0,
        (v_t + (v * v_x + g * h_x * dde.backend.backend.sin(alpha + np.pi/2))) - (g * h * dde.backend.backend.sin(alpha) * (1 - alpha_x) - C_f * v * v / h)
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
        num_boundary=100, 
        num_initial=160,
    )
    
    net = dde.nn.FNN(
        [2] + [20] * 3 + [20] * 3 + [20] * 3 + [2],  # Neural network structure
        'tanh',                # Activation function
        'Glorot normal',       # Weight initialization scheme
    )
    model = dde.Model(data, net)
    
    model.compile('adam', lr=1.0e-3)
    model.train(iterations=50000)
    model.compile('L-BFGS')
    losshistory, train_state = model.train()
    
    X = geomtime.random_points(100000)
    
    f = model.predict(X, operator=pde)
    err_eq = np.absolute(f)
    err = np.mean(err_eq)
    print("Mean residual: %.3e" % (err))
    
    dde.saveplot(losshistory, train_state, issave=True, isplot=True)
    
    samples = load_data(1000)
    y_true = samples[:, 2:4]
    y_pred = model.predict(samples[:, 0:2])
    print('L2 relative error:', dde.metrics.l2_relative_error(samples[:, 2:4], y_pred))
    exit()
    
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
    