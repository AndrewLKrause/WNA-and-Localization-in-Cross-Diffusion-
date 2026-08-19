parameters={} # Dictionary to store fixed parameters together with their values

unevaluatedparameters = [] # List to store the parameters that are to be kept symbolically

par1 = '' # First parameter in the expansion
par2 = '' # Second parameter in the expansion

translate = True # Boolean variable to state whether you want the system to be translated to the origin or not

equilibrium = [] # List of the strings representing the coordinates of the homogeneous steady state

var = [] # List of strings with the variables of the system

diffmatrix = [] List of lists to represent the diffusion matrix

kinetics = [] # List of functions in string form to repreesent the kinetics of the system

tol = 1e-7 # Tolerance value for the determinant of the jacobian matrix in order todetermine an invertible submatri to find the kernel of the noninvertible jacobian matrix

parameter_functions = {} # Dictionary that sets parameters as functions of other parameters rather than simply being evaluated

phiunit = 'n' # Boolean variable to state whether one wants the vector spanning the kernel of the jacobian matrix to be a unit vector or not

parameters_on_axes = [] # List with two strings representing the parameters that are used to compute bifurcation curves

names_of_parameters = [] # List to determine the names that are needed to be shown in the axes of the bifurcation diagram, if different from the names of the parameters used in the set up

intervalx = [] # List of two numbers in increasing order that set up the limits of the first parameter in the bifurcation diagram

intervaly = [] # List of two numbers in increasing order that set up the limits of the second parameter in the bifurcation diagram

lines_to_search = {} # Dictionary to set up the values of parameters to fix in order to find bifurcation curves. These values are only used to obtain initial conditions to find the bifurcation curves. These values do not remain fixed.