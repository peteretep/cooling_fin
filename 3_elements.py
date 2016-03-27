import numpy as np

grid_size = 5
l = 1.5
nsquared = 25
n = np.sqrt(nsquared)
Tb = 250
Tinf = 20

dx = l/grid_size

def finite_volume_method():
    a_matrix = []
    b_matrix =[]
    for x in range(grid_size):
        list = [0]*grid_size
        if x == 0:
            list[x] = 3/dx + nsquared * dx
            list[x+1] = -1/dx
            b = 2*Tb/dx + nsquared * Tinf * dx
        elif x == grid_size - 1:
            list[x-1] = -1/dx
            list[x] = 1/dx + (nsquared * dx)
            b = nsquared * Tinf * dx
        else:
            list[x-1] = -1/dx
            list[x] = 2/dx + nsquared * dx
            list[x+1] = -1/dx
            b = nsquared * Tinf * dx
        a_matrix.append(list)
        b_matrix.append(b)
    a = np.array(a_matrix)
    b = np.array(b_matrix)
    temps = np.linalg.solve(a, b)
    return temps


def analytical_solution():
    analytical_solutions = []
    for x in range(grid_size):
        element_size = l/grid_size
        r = element_size/2 + element_size*(x)
        T = Tinf + (np.cosh(n*(l - r))*(Tb - Tinf) / np.cosh(n*l))
        analytical_solutions.append(T)
    return np.array(analytical_solutions)

def discretisation_error():
    c=[]
    for i in range(len(analytical_solution())):
        c.append(float(analytical_solution()[i] - finite_volume_method()[i]) * 100 /analytical_solution()[i])
    return c

print('Analytical Solution:')
print analytical_solution()

print('Finite Volume Solution:')
print finite_volume_method()

print('Percentage Difference (discretisation error)')
print discretisation_error()

