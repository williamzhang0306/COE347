import numpy as np
from scipy.interpolate import interp1d

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
    - x (np.array): an array of xi values between [x0, x_N+1] inclusive of boundaries
    - u (np.array): the approximation of u(xi) for each xi value in x
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

    # put boundary conditions in data
    x = np.concatenate(([x0], interior_x, [xf]), axis = 0)
    u = np.concatenate(([u0], interior_u, [uf]), axis = 0)

    return x,u

def problem_5_solver(x0, uprime0, xf, uf, N, f):
    '''
    Solves the boundary value problem:
    u''(x) = f(x)
    u'(x0) = u_prime_0
    u(xf) = u_f

    Method: Use a second order one sided approximation of u'(x0)
    and apply a second order central differene to find all internior u(x)

    Parameters
    - x0, the left boundary
    - u_prime_0: the boundary condition u'(x0)
    - xf, the right boundary
    - uf: the boundary condition u(xf)
    - N: the number of interior points
    - f: u''(x) = f(x)
    '''
    # working interval - [x0, xf) left boundary and interior points, execlude left boundary for now
    h = (xf - x0)/(N+1)
    x_values = np.linspace(x0, xf - h, N+1)

    # construct tridiagonal (N+1 x N+1) matrix
    main_diagonal = -2 * np.ones(N+1)
    above_below_diagonal = np.ones(N)
    A = np.diag(main_diagonal) + np.diag(above_below_diagonal, 1) + np.diag(above_below_diagonal, -1)

    # put one sided dfd approximation of u(0) in the top row
    A[0][0] = -3/2
    A[0][1] = 2
    A[0][2] = -1/2

    #lhs
    b = (h**2) * np.array([f(xi) for xi in x_values])
    b[0] = uprime0 * h
    b[-1] -= uf
    
    u = np.linalg.solve(A,b)

    # add last boundary
    x_values = np.append(x_values,xf)
    u = np.append(u,uf)

    return x_values,u

def problem_6_solver(x0, uprime0, xf, uf, N, f):
    '''
    Solves the boundary value problem:
    u''(x) = f(x)
    u'(x0) = u_prime_0
    u(xf) = u_f

    Method: Use a first order one sided approximation of u'(x0)
    and apply a second order central differene to find all internior u(x)

    Parameters
    - x0, the left boundary
    - u_prime_0: the boundary condition u'(x0)
    - xf, the right boundary
    - uf: the boundary condition u(xf)
    - N: the number of interior points
    - f: u''(x) = f(x)
    '''
    # working interval - [x0, xf) left boundary and interior points, execlude left boundary for now
    h = (xf - x0)/(N+1)
    x_values = np.linspace(x0, xf - h, N+1)

    # construct tridiagonal (N+1 x N+1) matrix
    main_diagonal = -2 * np.ones(N+1)
    above_below_diagonal = np.ones(N)
    A = np.diag(main_diagonal) + np.diag(above_below_diagonal, 1) + np.diag(above_below_diagonal, -1)

    # put first order one sided dfd approximation of u(0) in the top row
    A[0][0] = -1
    A[0][1] = 1

    #lhs
    b = (h**2) * np.array([f(xi) for xi in x_values])
    b[0] = uprime0 * h
    b[-1] -= uf
    
    u = np.linalg.solve(A,b)

    # add last boundary
    x_values = np.append(x_values,xf)
    u = np.append(u,uf)

    return x_values,u

def error_1(x_list:np.array, u_list:np.array, exact_x_data:np.array, exact_u_data:np.array):
    '''
    Calculates the error denoted 'E', defined by problem 3

    Parameters:
    - xi (list or similar):
    - ui 
    - exact_xi
    - exact_ui

    Returns
    - Error (float)
    '''

    # interpolated function
    exact_u_interpolation = interp1d(exact_x_data,exact_u_data, kind='quadratic')

    error = sum([(u_xi - exact_u_interpolation(x_i))**2 for x_i, u_xi in zip(x_list,u_list)])

    return error**.5

def error_2(x_list:np.array, u_list:np.array, exact_x:np.array, exact_u:np.array, N):
    return (1/N) * error_1(x_list, u_list, exact_x, exact_u)



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
