import numpy as np

def ratio_test(A, b, x, p):
    feasible_directions = []
    for i, (a_i, b_i) in enumerate(zip(A, b)):
        feasible = np.less(np.matmul(a_i, p), 0)
        feasible_directions.append((i, np.matmul(a_i, p)))
    
        if feasible:
            ratio = (np.matmul(a_i, x) - b_i)/np.matmul(-a_i, p)
            print('a_{} ratio: ({} - {})/{} = '.format(i, np.matmul(a_i, x), b_i, np.matmul(-a_i, p)), ratio)
    
    [print('a_{}^T p: {}'.format(*d)) for d in feasible_directions]