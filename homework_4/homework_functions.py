import numpy as np

def example_f(x):
    '''f(x) = u''(x) for all problems in homework 4'''
    return (4 * np.pi)**2 * np.cos(4*np.pi * x)


def centered_diff(boundary1:tuple , boundary2:tuple,N:int,f):
    '''
    Uses a second order centered difference approximation
    of the second derivative to compute function values ( u(x) ) within an interval.

    Paramters:
    - boundary1 (tuple): (x0, u(x0)) - the boundary condition at x = x0
    - boundary2 (tuple): (x_N+1, u(x_N+1)) - the boundary condition at x = x_N+1
    - N (int): the number of interior points in the interval [x0,x_N+1]
    - f (func): a function f(x) = u''(x)

    Returns:
    - interior_x (np.array): an array of interior xi values between (x0, x_N+1)
    - interior_u (np.array): the approximation of u(xi) for each xi value in interior_x
    '''

    # get interval
    x0, u0 = boundary1
    xf, uf = boundary2
    h = (xf - x0)/(N+1)
    interior_x = np.linspace(x0 + h , xf - h, N)

    # construct tridiagonal matrix
    main_diagonal = -2 * np.ones(N)
    above_below_diagonal = np.ones(N-1)
    A = np.diag(main_diagonal) + np.diag(above_below_diagonal, 1) + np.diag(above_below_diagonal, -1)

    # LHS - h^2f + boundary conditions
    b = (h**2) * np.array([f(xi) for xi in interior_x])
    b[0] -= u0
    b[-1] -= uf

    # solve for u at interior points
    interior_u = np.linalg.solve(A, b)

    return interior_x, interior_u


def dat_reader(file_path):
    '''
    Reads '.dat' files into numpy arrays

    Paramters:
    - file_path (str): .dat file to open

    Returns:
    - first_column (np.array): the first column of data in the .dat file
    - second_column (np.array): the second column of data in the .dat file
    '''
    # Read the values from the file
    data = np.loadtxt(file_path)

    # Separate the values into two arrays
    first_column = data[:, 0]
    second_column = data[:, 1]

    return first_column, second_column
