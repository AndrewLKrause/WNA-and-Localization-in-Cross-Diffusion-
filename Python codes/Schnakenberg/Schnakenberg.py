parameters = {'d': 1}

unevaluatedparameters = []

par1 = 'b'
par2 = 'p'

translate = True

equilibrium = ['b', '1/b']

var = ['u',
       'v']

diffmatrix = [[1, 'p*u/(1 + u**2)'],
              [0, 'd']]

kinetics = ['- u + u**2*v',
            'b - u**2*v']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

parameters_on_axes = ['b', 'p']

names_of_parameters = ['\\lambda', 'p']

intervalx = [0.4, 1.3]

intervaly = [0, 2.5]

lines_to_search = {'p':2, 'b':1.05}