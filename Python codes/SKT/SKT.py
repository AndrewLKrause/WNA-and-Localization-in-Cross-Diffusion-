parameters={'r1': 1.0,
            'r2': 2.0,
            'a1': 0.9,
            'b1': 0.6,
            'b2': 0.2,
            'd1': 1.0,
            'd2': 1.0,
            'd11': 0.0,
            'd12': 0.0,
            'd22': 0.0}

unevaluatedparameters = []

par1 = 'p'
par2 = 'q'

translate = True

equilibrium = ['(p - b1)/(p*a1 - b1*b2)', '(a1 - b2)/(p*a1 - b1*b2)']

var = ['u',
       'v']

diffmatrix = [['d1 + d11*u + d12*v', 'd12*u'],
              ['q*v', 'd2 + q*u + d22*v']]

kinetics = ['r1*u*(1 - a1*u - b1*v)',
            'r2*v*(1 - b2*u - p*v)']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

parameters_on_axes = ['p', 'q']

names_of_parameters = ['a_2', 'd_{21}']

intervalx = [0.5, 1.0]

intervaly = [0, 200]

lines_to_search = {'p': 0.9}