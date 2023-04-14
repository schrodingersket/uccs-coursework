import numpy as np

def gll_newton_raphson(N, tol=10e-8):
    """
    Computes N Gauss-Legendre-Lobatto nodes (intervals) on the interval [-1, 1] by using the Newton-Raphson method to 
    compute zeros of (1 - x^2)P_N'(x), where P_N' denotes the first derivative of the N-th degree Legendre polynomial.
    
    To use the Newton-Raphson method, we need d/dt [(1 - x^2)P_N'(x)]. 
    
    See https://doi.org/10.1007/978-88-470-5522-3_10
    """

    # Initial guess
    #
    x = np.polynomial.chebyshev.chebpts2(N+1)

    # Legendre polynomials
    #
    L_N = np.array((
        np.polynomial.legendre.Legendre.basis(N),
        np.polynomial.legendre.Legendre.basis(N-1),
    ))

    # Generate Legendre polynomials
    #
    x_prev = x
    for i in range(100):
        x = x + (L_N[-1](x) - x * L_N[0](x)) / ((N + 1) * L_N[0](x))

        if np.linalg.norm(x - x_prev, np.inf) < tol:
            break
    
    # Compute weights
    # 
    w = 2.0 / ((N*(N + 1)) * (L_N[0](x)**2))

    return x, w

def gll(N, *args):
    """
    Computes N Gauss-Legendre-Lobatto nodes (intervals) on the interval [-1, 1] by using Numpy's built-in libraries for
    finding roots of differentiated Legendre polynomials.
    """
    # Get n-th Legendre polynomial basis function
    #
    L_N = np.polynomial.legendre.Legendre.basis(N)

    # Find zeros of derivative of Legendre polynomials
    #
    x = np.polynomial.Legendre.basis(N).deriv().roots()
    x = np.insert(x, 0, -1.0)
    x = np.append(x, 1.0)

    # Compute weights
    # 
    w = 2.0 / ((N*(N + 1)) * (L_N(x)**2))

    return x, w

if __name__ == '__main__':
    import time
    np.set_printoptions(suppress = True)
    print('Gauss-Legendre-Lobatto nodes/weights for 4th-degree Legendre polynomial:')

    x4, w4 = gll(4)
    print('x = {}'.format(x4))
    print('w = {}'.format(w4))

    # Compute GLL via Newton-Raphson
    #
    start = time.time()
    x_nr, w_nr = gll_newton_raphson(100, 10e-8)
    end = time.time()
    print('Newton-Raphson took {} seconds'.format(end - start))
    start = time.time()

    # Compute GLL via Numpy
    #
    x_np, w_np = gll(100)
    end = time.time()
    print('Numpy differentiation took {} seconds'.format(end - start))

    # Assert equivalence
    #
    assert(np.all(np.isclose(x_nr, x_np)))
    assert(np.all(np.isclose(w_nr, w_np)))