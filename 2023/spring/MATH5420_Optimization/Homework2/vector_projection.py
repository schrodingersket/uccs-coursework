import numpy as np

def proj(u, v):
    return ((np.dot(u, v)/(np.linalg.norm(v)**2)), v)