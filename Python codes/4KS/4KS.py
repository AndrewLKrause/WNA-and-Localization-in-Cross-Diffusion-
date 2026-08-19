parameters={'a': 0.9,
            'D': 1.0,
            'r0': 0.1,
            'p': 0.2,
            'q': 0.5}

unevaluatedparameters = []

par1 = 'b'
par2 = 'c'

translate = True

equilibrium = [1, '(b + p)/(a*(b + p) + p*q)', 1, 'q/(a*(b + p) + p*q)']

var = ['u',
       'v',
       'w',
       'z']

diffmatrix = [[1, '- c*u/(1 + u**2)', 0, 0],
              [0, 'D', 0, 0],
              [0, 0, 'D', 0],
              [0, 0, 0, 'D']]

kinetics = ['r0*u*(1 - u)',
            'u - a*v - q*v*w + b*z*u**2',
            '- q*v*w + (p + b*u**2)*z + u*(1 - w)',
            'q*v*w - (p + b*u**2)*z']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

parameters_on_axes = ['b', 'c']

names_of_parameters = []

intervalx = [0, 2.0]

intervaly = [2.5, 4.5]

lines_to_search = {'b': 0.1}