parameters={'k': 0.2,
            'a': 20.0,
            'D': 0.2}

par1 = 'r'
par2 = 'p'

unevaluatedparameters = {}

translate = True

equilibrium = ['1', '1/a']

var = ['u',
       'v']

diffmatrix = [['1', '- p*u/(1 + u**2)'],
              [0, 'D']]

kinetics = ['r*u*(1 - u)/(k + v)',
            'u - a*v']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

parameters_on_axes = ['r', 'p']

names_of_parameters = []

intervalx = [0, 5]

intervaly = [50, 100]

lines_to_search = {'r': 1.0}