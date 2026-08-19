from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf')

import os
import shutil
from sys import exit
from sympy import *
from sympy.printing.mathematica import mathematica_code
init_printing()
from sympy.solvers import solve
from mpmath import findroot
import math

# Setting up the variables for the expansion. This includes shifting system to the origin and checking that the steady state provided are correct

# modelname = '4KS'
# modelname = 'HdW2'
# modelname = 'Keller-Segel - C2'
# modelname = 'SKT'
# modelname = 'Schnakenberg'
modelname = input('Enter the name of the system you would like to analyze: ')

# Checking for directory

if not os.path.isdir(modelname):
    print('The directory of that system was not found. Create it first and place ' + 
          'the data file inside named in the same way.')
    exit()
else:
    os.chdir(os.getcwd() + '\\' + modelname)

try:
    exec(open(modelname + '.py').read())
except:
    print('The file ' + modelname + '.py could not be run')
    exit()
    
nvar = len(var)

# Checking for the file functions.py

try:
    exec(open(os.path.dirname(os.path.realpath(__file__)) + '\\functions.py').read())
except:
    print('File functions.py is not in the same folder as the script you are running')
    exit()

file = open('List of Variables.txt', 'w')

# Checking for variables of the system

try:
    var
except:
    print('Variables were not provided')
    exit()

for varnum in range(nvar):
    try:
        exec(var[varnum] + ' = symbols(var[varnum], real=True)')
        var[varnum] = eval(var[varnum])
    except:
        print('The script could not define ' + (var[varnum]) + ' as a variable')
        exit()
    file.write(mathematica_code(var[varnum]) + '\n')

# Checking for parameters

try:
    parameters
except:
    print('Parameters were not provided. The script will assume that there are no parameters.')
    parameters=[]

npar = len(parameters) + len(unevaluatedparameters)

if npar>0:
    newparameters = dict()
    for key in parameters.keys():
        try:
            exec(key + ' = symbols(key, real = True)')
            newparameters[eval(key)] = parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters = newparameters

for parnum in range(len(unevaluatedparameters)):
    exec(unevaluatedparameters[parnum] + ' = symbols(unevaluatedparameters[parnum], real=True)')
    unevaluatedparameters[parnum] = eval(unevaluatedparameters[parnum])
    
# Checking for diffusion matrix

try:
    diffmatrix
except:
    print('The diffusion matrix was not provided correctly.')
    exit()
    
exec(par1 + ' = symbols(par1, real = True)')
exec(par2 + ' = symbols(par2, real = True)')

par1 = eval(par1)
par2 = eval(par2)

# Checking for equilibrium

for varnum in range(nvar):
    try:
        equilibrium[varnum] = eval(equilibrium[varnum])
    except:
        pass

try:
    for row in range(nvar):
        for col in range(nvar):
            diffmatrix[row][col] = eval(str(diffmatrix[row][col]).replace('^', '**'))
except:
    print('The diffusion matrix is not a function of the parameters of the system')
    exit()
        
diffmatrix = Matrix(diffmatrix)

# Checking for kinetics

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()
    
for functionnumber in range(nvar):
    try:
        kinetics[functionnumber] = eval(str(kinetics[functionnumber]).replace('^', '**'))
    except:
        print('The expression ' + kinetics[functionnumber] + ' is not a function of the parameters of your system')
        exit()

kinetics = Matrix(kinetics)

if modelname!='4SKT' and modelname!='Andrew':
    if not check_equilibrium(equilibrium, kinetics):
        print('The coordinates of the equilibrium do not solve the trivial system.')
        exit()

if equilibrium!=[0]*nvar:
    equilibrium, kinetics, diffmatrix = origin_translation(equilibrium, kinetics, diffmatrix)
    
eq_subs = dict(zip(var, equilibrium))
    
aNF = []
bNF = []

for counter in range(5):
    aNF.append(symbols(str(par1) + 'NF' + str(counter), real = True))
    bNF.append(symbols(str(par2) + 'NF' + str(counter), real = True))

extraparvals = {}

kinetics = parameter_translation(kinetics, par1, aNF[0], par2, bNF[0])
diffmatrix = parameter_translationdiff(diffmatrix, par1, aNF[0], par2, bNF[0])

expandingparvals = {par1: 0, par2: 0}

file = open('Kinetics.txt', 'w')
for functionnumber in range(nvar):
    file.write(mathematica_code(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat = kinetics.jacobian(var).subs(expandingparvals)
jacobianmat = jacobianmat.subs(eq_subs)

# Setting up variables like the wavenumber and derivatives of the vector field to generate the vectors defined via the asymptotic expansion

muNF = symbols(r'muNF', real = true)

for parnum1 in range(5):
    for parnum2 in range(5):
        exec('f' + str(parnum1) + str(parnum2) + '= Mul(Pow(Mul(math.factorial(parnum1), ' + 
             'math.factorial(parnum2)), -1), diff(kinetics, par1, parnum1, par2, parnum2))' +
             '.subs(par1, 0).subs(par2, 0)')
        
for parnum1 in range(5):
    for parnum2 in range(5):
        exec('diffmatrix' + str(parnum1) + str(parnum2) + '= Mul(Pow(Mul(math.factorial(parnum1), ' + 
             'math.factorial(parnum2)), - 1), diff(diffmatrix, par1, parnum1, par2, parnum2))' +
             '.subs(par1, 0).subs(par2, 0)')

for counter in range(5):
    exec(f'coefmat{counter} = matrix("coef{counter}mat", nvar' + \
         ', Add(jacobianmat, Mul(- Pow(counter, 2), muNF, diffmatrix00.subs(eq_subs))))')

# Turing bifurcation condition

TC1 = coefmat1.dummy.det()

for row in range(nvar):
    for col in range(nvar):
        TC1 = TC1.subs(coefmat1.dummy[row, col], coefmat1.actualcoord[row, col])
        
TC2 = diff(TC1, muNF)

eigval = symbols('eigval', real=True)

auxmat = matrix('auxmat', nvar, Add(jacobianmat, Mul(- 1, muNF, diffmatrix), Mul(- 1, eigval, eye(nvar))))

auxterm = auxmat.dummy.det()

for row in range(nvar):
    for col in range(nvar):
        auxterm = auxterm.subs(auxmat.dummy[row, col], auxmat.actualcoord[row, col])

nonzero = diff(auxterm, eigval).subs(eigval, 0)
        
file = open('Determinant of the Jacobian matrix.txt', 'w')
file.write(mathematica_code(TC1))
file.close()
file = open('Derivative of the Determinant.txt','w')
file.write(mathematica_code(TC2))
file.close()
file = open("Non-zero variable.txt", 'w')
file.write(mathematica_code(nonzero))
file.close()

negativeRHS = Vector('negativeRHS')

# Computing derivatives of the vector field up to order 5

for parnum1 in range(5):
    for parnum2 in range(5):
        if parnum1 + parnum2<=6:
            exec(f'firstorderderivatives{parnum1}{parnum2} = list()')
            exec(f'firstorderderivativesdiff{parnum1}{parnum2} = list()')
            if parnum1 + parnum2<=5:
                exec(f'secondorderderivatives{parnum1}{parnum2} = list()')
                exec(f'secondorderderivativesdiff{parnum1}{parnum2} = list()')
                if parnum1 + parnum2<=4:
                    exec(f'thirdorderderivatives{parnum1}{parnum2} = list()')
                    exec(f'thirdorderderivativesdiff{parnum1}{parnum2} = list()')
                    if parnum1 + parnum2<=3:
                        exec(f'fourthorderderivatives{parnum1}{parnum2} = list()')
                        exec(f'fourthorderderivativesdiff{parnum1}{parnum2} = list()')
                        if parnum1 + parnum2<=2:
                            exec(f'fifthorderderivatives{parnum1}{parnum2} = list()')
        for counter1 in range(nvar):
            if parnum1 + parnum2<=6:
                exec(f'auxderc1 = diff(f{parnum1}{parnum2}, var[{counter1}])')
                exec(f'firstorderderivatives{parnum1}{parnum2}.append(auxderc1.subs(eq_subs))')
                exec(f'auxderc1diff = diff(diffmatrix{parnum1}{parnum2}, var[{counter1}])')
                exec(f'firstorderderivativesdiff{parnum1}{parnum2}.append(auxderc1diff.subs(eq_subs))')
                if parnum1 + parnum2<=5:
                    exec(f'secondorderderivatives{parnum1}{parnum2}.append(list())')
                    exec(f'secondorderderivativesdiff{parnum1}{parnum2}.append(list())')
                    if parnum1 + parnum2<=4:
                        exec(f'thirdorderderivatives{parnum1}{parnum2}.append(list())')
                        exec(f'thirdorderderivativesdiff{parnum1}{parnum2}.append(list())')
                        if parnum1 + parnum2<=3:
                            exec(f'fourthorderderivatives{parnum1}{parnum2}.append(list())')
                            exec(f'fourthorderderivativesdiff{parnum1}{parnum2}.append(list())')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}.append(list())')
            for counter2 in range(nvar):
                if parnum1 + parnum2<=5:
                    exec(f'auxderc2 = diff(auxderc1, var[{counter2}])')
                    exec(f'secondorderderivatives{parnum1}{parnum2}[{counter1}].append(auxderc2.subs(eq_subs))')
                    exec(f'auxderc2diff = diff(auxderc1diff, var[{counter2}])')
                    exec(f'secondorderderivativesdiff{parnum1}{parnum2}[{counter1}].append(auxderc2diff.subs(eq_subs))')
                    if parnum1 + parnum2<=4:
                        exec(f'thirdorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                        exec(f'thirdorderderivativesdiff{parnum1}{parnum2}[{counter1}].append(list())')
                        if parnum1 + parnum2<=3:
                            exec(f'fourthorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                            exec(f'fourthorderderivativesdiff{parnum1}{parnum2}[{counter1}].append(list())')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                for counter3 in range(nvar):
                    if parnum1 + parnum2<=4:
                        exec(f'auxderc3 = diff(auxderc2, var[{counter3}])')
                        exec(f'thirdorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(auxderc3.subs(eq_subs))')
                        exec(f'auxderc3diff = diff(auxderc2diff, var[{counter3}])')
                        exec(f'thirdorderderivativesdiff{parnum1}{parnum2}[{counter1}][{counter2}].append(auxderc3diff.subs(eq_subs))')
                        if parnum1 + parnum2<=3:
                            exec(f'fourthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                            exec(f'fourthorderderivativesdiff{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                    for counter4 in range(nvar):
                        if parnum1 + parnum2<=3:
                            exec(f'auxderc4 = diff(auxderc3, var[{counter4}])')
                            exec(f'fourthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(auxderc4.subs(eq_subs))')
                            exec(f'auxderc4diff = diff(auxderc3diff, var[{counter4}])')
                            exec(f'fourthorderderivativesdiff{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(auxderc4diff.subs(eq_subs))')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(list())')
                        for counter5 in range(nvar):
                            if parnum1 + parnum2<=2:
                                exec(f'auxderc5 = diff(auxderc4, var[{counter5}])')
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}].append(auxderc5.subs(eq_subs))')

for counter in range(1, 5):
    exec(f'auxderc{counter} = 0')
    exec(f'auxderc{counter}diff = 0')
    
# Symbolic definition of the vectors present in the expansion

phiNF = Vector('phiNF')
psiNF = Vector('psiNF')

# Order 2 vectors

for counter1 in range(0, 3, 2):
    exec(f'W{counter1}2NF = Vector("W{counter1}2NF")')
    
# Order 3 vectors
    
for counter1 in range(1, 4, 2):
    if counter1!=1:
        exec(f'W{counter1}3NF = Vector("W{counter1}3NF")')
    else:
        for counter2 in range(1, 4):
            if counter2==1:
                exec(f'W{counter1}3NF = Vector("W{counter1}3NF")')
            else:
                exec(f'W{counter1}{counter2}3NF = Vector("W{counter1}{counter2}3NF")')

# Order 4 vectors
        
for counter1 in range(0, 3, 2):
    if counter1<3:
        for counter2 in range(1, 4):
            if counter2==1:
                exec(f'W{counter1}4NF = Vector("W{counter1}4NF")')
            else:
                exec(f'W{counter1}{counter2}4NF = Vector("W{counter1}{counter2}4NF")')
    else:
        exec(f'W{counter1}4NF = Vector("W{counter1}4NF")')

getout = 0

# Determination of a row and a column that can be removed from noninvertible jacobian matrix to obtain an invertible submatrix

for row in range(nvar):
    for col in range(nvar):
        submatrixrows = list(range(nvar))
        submatrixcols = list(range(nvar))
        submatrixrows.remove(row)
        submatrixcols.remove(col)
        invertiblesubmatrix = coefmat1.actualcoord.extract(submatrixrows, submatrixcols)
        submatrixeval = invertiblesubmatrix
        submatrixeval = submatrixeval.subs(parameters)
        for varnum in range(nvar):
            submatrixeval = submatrixeval.subs(var[varnum], equilibrium[varnum])
        try:
            if abs(N(submatrixeval.det())) > tol:
                phiNF.actualcoord[col] = 1
                criticalrow = row
                criticalcol = col
                getout = 1
                break
        except:
            phiNF.actualcoord[col] = 1
            criticalrow = row
            criticalcol = col
            getout = 1
            break
    if getout==1:
        break
    
# Determination of the kernel of the Jacobian matrix and its transpose
    
coefsubmatrix = matrix('dummysubmatrix', nvar - 1, invertiblesubmatrix)

phiNF = kernel_determination(phiNF, coefmat1, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

if phiunit=='y':
    phiNF.actualcoord = Mul(Pow(sqrt(phiNF.actualcoord.dot(phiNF.actualcoord)), - 1), phiNF.actualcoord)

Tcoefsubmatrix = matrix('Tdummysubmatrix', nvar - 1, transpose(invertiblesubmatrix))
Tcoefmat1 = matrix('Tcoefmat', nvar, transpose(coefmat1.actualcoord))

psiNF = kernel_determination(psiNF, Tcoefmat1, criticalrow, Tcoefsubmatrix, submatrixcols, submatrixrows)

for row in range(nvar):
    for col in range(nvar):
        psiNF.actualcoord = psiNF.actualcoord.subs(coefmat1.dummy[row, col], coefmat1.actualcoord[row, col])
        if row<nvar-1 and col<nvar-1:
            psiNF.actualcoord = psiNF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                       coefsubmatrix.actualcoord[row, col])

phiNF_eval = evaluation_dict(phiNF)
psiNF_eval = evaluation_dict(psiNF)

print('First order ready')

# In what comes below, the right-hand side of the equation provided in the paper is defined in the variable negativeRHS.actualcoord and then the corresponding system is solved.
# When the system has a non-invertible matrix of coefficients, an equation is set up by using the Fredholm alternative. The corresponding equations are solved in order to find
# the right coefficients of the expansion

# Specifically, the variable defined below represents the second order Taylor expansion of f_{0, 0} at 0 in the directions phi, phi. Once again, this logic is repeated throughout
# what comes below

DS_phiphi_00 = second_order(0, 0, phiNF, phiNF)

negativeRHS.actualcoord = DS_phiphi_00

W02NF = linearsolver(W02NF, negativeRHS, coefmat0)

SSD_phi_00 = first_order_diff(0, 0, phiNF)

negativeRHS.actualcoord = Add(DS_phiphi_00, Mul(- 2, muNF, SSD_phi_00, phiNF.dummy))

W22NF = linearsolver(W22NF, negativeRHS, coefmat2)

for counter in range(0, 3, 2):
    exec(f'W{counter}2NF_eval = evaluation_dict(W{counter}2NF)')

print('Second order ready')

SS_phi_10 = first_order(1, 0, phiNF)
SS_phi_01 = first_order(0, 1, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[2], SS_phi_10), Mul(bNF[2], SS_phi_01),
                              Mul(- muNF, Add(Mul(aNF[2], diffmatrix10),
                                              Mul(bNF[2], diffmatrix01)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

upsilon1 = psiNF.dummy.dot(negativeRHS.actualcoord)

W13NF = critical_linearsolver(W13NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

DS_phiW02_00 = second_order(0, 0, phiNF, W02NF)
DS_phiW22_00 = second_order(0, 0, phiNF, W22NF)
TS_phiphiphi_00 = third_order(0, 0, phiNF, phiNF, phiNF)

SSD_W02_00 = first_order_diff(0, 0, W02NF)
SSD_W22_00 = first_order_diff(0, 0, W22NF)
DSD_phiphi_00 = second_order_diff(0, 0, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(2, Add(Mul(2, DS_phiW02_00), DS_phiW22_00)),
                              Mul(3, TS_phiphiphi_00),
                              Mul(- 2, muNF, SSD_phi_00, W22NF.dummy),
                              Mul(- muNF, Add(Mul(2, SSD_W02_00), - SSD_W22_00,
                                              DSD_phiphi_00), phiNF.dummy)).subs(extraparvals)

# Codimension-two condition

upsilon2 = psiNF.dummy.dot(negativeRHS.actualcoord)

file = open('Third-order coefficient.txt', 'w')
file.write(mathematica_code(upsilon2))
file.close()
    
W123NF = critical_linearsolver(W123NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

negativeRHS.actualcoord = Mul(2, sqrt(muNF), diffmatrix00, phiNF.dummy).subs(eq_subs)

W133NF = critical_linearsolver(W133NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW22_00), TS_phiphiphi_00,
                              Mul(- 6, muNF, SSD_phi_00, W22NF.dummy),
                              Mul(- 3, muNF, Add(SSD_W22_00, DSD_phiphi_00),
                                  phiNF.dummy)).subs(extraparvals)

W33NF = linearsolver(W33NF, negativeRHS, coefmat3)

for counter in range(1, 4, 2):
    exec(f'W{counter}3NF_eval = evaluation_dict(W{counter}3NF)')
    
for counter in range(2, 4):
    exec(f'W1{counter}3NF_eval = evaluation_dict(W1{counter}3NF)')

print('Third order ready')

SS_W02_10 = first_order(1, 0, W02NF)
SS_W02_01 = first_order(0, 1, W02NF)
DS_phiW13_00 = second_order(0, 0, phiNF, W13NF)
DS_phiphi_10 = second_order(1, 0, phiNF, phiNF)
DS_phiphi_01 = second_order(0, 1, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[2], SS_W02_10), Mul(bNF[2], SS_W02_01),
                              Mul(2, DS_phiW13_00), Mul(aNF[2], DS_phiphi_10),
                              Mul(bNF[2], DS_phiphi_01)).subs(extraparvals)

W04NF = linearsolver(W04NF, negativeRHS, coefmat0)

DS_phiW123_00 = second_order(0, 0, phiNF, W123NF)
DS_W02W02_00 = second_order(0, 0, W02NF, W02NF)
DS_W22W22_00 = second_order(0, 0, W22NF, W22NF)
TS_phiphiW02_00 = third_order(0, 0, phiNF, phiNF, W02NF)
TS_phiphiW22_00 = third_order(0, 0, phiNF, phiNF, W22NF)
Q4S_phiphiphiphi_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW123_00), Mul(2, DS_W02W02_00),
                              DS_W22W22_00, Mul(3, Add(Mul(2, TS_phiphiW02_00), TS_phiphiW22_00)),
                              Mul(3, Q4S_phiphiphiphi_00))

W024NF = linearsolver(W024NF, negativeRHS, coefmat0)

DS_phiW133_00 = second_order(0, 0, phiNF, W133NF)

negativeRHS.actualcoord = Mul(2, DS_phiW133_00)

W034NF = linearsolver(W034NF, negativeRHS, coefmat0)

SS_W22_10 = first_order(1, 0, W22NF)
SS_W22_01 = first_order(0, 1, W22NF)

SSD_phi_10 = first_order_diff(1, 0, phiNF)
SSD_phi_01 = first_order_diff(0, 1, phiNF)
SSD_W13_00 = first_order_diff(0, 0, W13NF)

negativeRHS.actualcoord = Add(Mul(aNF[2], SS_W22_10), Mul(bNF[2], SS_W22_01), Mul(2, DS_phiW13_00),
                              Mul(aNF[2], DS_phiphi_10), Mul(bNF[2], DS_phiphi_01),
                              Mul(- 2, muNF, SSD_phi_00, W13NF.dummy),
                              Mul(- 4, muNF, Add(Mul(aNF[2], diffmatrix10),
                                                 Mul(bNF[2], diffmatrix01)),
                                  W22NF.dummy),
                              Mul(- 2, muNF, Add(SSD_W13_00,
                                                 Mul(aNF[2], SSD_phi_10),
                                                 Mul(bNF[2], SSD_phi_01)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

W24NF = linearsolver(W24NF, negativeRHS, coefmat2)

DS_phiW33_00 = second_order(0, 0, phiNF, W33NF)
DS_W02W22_00 = second_order(0, 0, W02NF, W22NF)

SSD_W123_00 = first_order_diff(0, 0, W123NF)
SSD_W33_00 = first_order_diff(0, 0, W33NF)
DSD_phiW02_00 = second_order_diff(0, 0, phiNF, W02NF)
TSD_phiphiphi_00 = third_order_diff(0, 0, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(2, Add(DS_phiW123_00, DS_phiW33_00)), Mul(4, DS_W02W22_00),
                              Mul(6, Add(TS_phiphiW02_00, TS_phiphiW22_00)), Mul(4, Q4S_phiphiphiphi_00),
                              Mul(- 2, muNF, SSD_phi_00, Add(W123NF.dummy,
                                                             Mul(3, W33NF.dummy))),
                              Mul(- 8, muNF, Add(SSD_W02_00, DSD_phiphi_00), W22NF.dummy),
                              Mul(- 2, muNF, Add(SSD_W123_00, - SSD_W33_00,
                                                 Mul(4, DSD_phiW02_00),
                                                 Mul(2, TSD_phiphiphi_00)),
                                  phiNF.dummy)).subs(extraparvals)

W224NF = linearsolver(W224NF, negativeRHS, coefmat2)

SSD_W133_00 = first_order_diff(0, 0, W133NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW133_00), Mul(- 2, muNF, SSD_phi_00, W133NF.dummy),
                              Mul(- 2, muNF, SSD_W133_00, phiNF.dummy),
                              Mul(8, sqrt(muNF), diffmatrix00, W22NF.dummy),
                              Mul(4, sqrt(muNF), SSD_phi_00, phiNF.dummy)).subs(eq_subs).subs(extraparvals)

W234NF = linearsolver(W234NF, negativeRHS, coefmat2)

for counter1 in range(0, 3, 2):
    exec(f'W{counter1}4NF_eval = evaluation_dict(W{counter1}4NF)')
    for counter2 in range(2, 4):
        exec(f'W{counter1}{counter2}4NF_eval = evaluation_dict(W{counter1}{counter2}4NF)')
    
print('Fourth order ready')

# Determination of the coefficients of the fifth-order amplitude equation and the corresponding vectors that solve the fifth order system

for counter in range(1, 4, 2):
    exec(f'W{counter}3NF_eval = evaluation_dict(W{counter}3NF)')

W04NF_eval = evaluation_dict(W04NF)
W24NF_eval = evaluation_dict(W24NF)
    
negativeRHS.actualcoord = Add(Mul(diffmatrix00, phiNF.dummy),
                              Mul(- 2, sqrt(muNF),
                                  diffmatrix00, W133NF.dummy)).subs(eq_subs)

alpha1 = psiNF.dummy.dot(negativeRHS.actualcoord)

SS_W133_10 = first_order(1, 0, W133NF)
SS_W133_01 = first_order(0, 1, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[2], SS_W133_10), Mul(bNF[2], SS_W133_01),
                              Mul(2, sqrt(muNF), diffmatrix00, W13NF.dummy),
                              Mul(- muNF, Add(Mul(aNF[2], diffmatrix10),
                                              Mul(bNF[2], diffmatrix01)),
                                  W133NF.dummy),
                              Mul(2, sqrt(muNF), Add(Mul(aNF[2], diffmatrix10),
                                                     Mul(bNF[2], diffmatrix01)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

alpha2 = psiNF.dummy.dot(negativeRHS.actualcoord)

DS_phiW034_00 = second_order(0, 0, phiNF, W034NF)
DS_phiW234_00 = second_order(0, 0, phiNF, W234NF)
DS_W02W133_00 = second_order(0, 0, W02NF, W133NF)
TS_phiphiW133_00 = third_order(0, 0, phiNF, phiNF, W133NF)

SSD_W034_00 = first_order_diff(0, 0, W034NF)
SSD_W234_00 = first_order_diff(0, 0, W234NF)

negativeRHS.actualcoord = Add(Mul(2, Add(DS_phiW034_00, DS_phiW234_00)),
                              Mul(4, DS_W02W133_00), Mul(6, TS_phiphiW133_00),
                              Mul(- 2, muNF, SSD_phi_00, W234NF.dummy),
                              Mul(4, sqrt(muNF), diffmatrix00, W123NF.dummy),
                              Mul(- 2, muNF, Add(SSD_W02_00, DSD_phiphi_00),
                                  W133NF.dummy),
                              Mul(2, sqrt(muNF), SSD_phi_00, Add(W02NF.dummy, Mul(3, W22NF.dummy))),
                              Mul(- muNF, Add(SSD_W034_00, - SSD_W234_00), phiNF.dummy),
                              Mul(2, sqrt(muNF), Add(Mul(3, SSD_W02_00),
                                                     - SSD_W22_00,
                                                     Mul(2, DSD_phiphi_00)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

alpha3 = psiNF.dummy.dot(negativeRHS.actualcoord)

DS_W22W133_00 = second_order(0, 0, W22NF, W133NF)

DSD_phiW133_00 = second_order_diff(0, 0, phiNF, W133NF)

negativeRHS.actualcoord =  Add(Mul(- 2, DS_phiW034_00), Mul(- 2, DS_W22W133_00),
                               Mul(- 3, TS_phiphiW133_00),
                               Mul(2, sqrt(muNF), diffmatrix00, W123NF.dummy),
                               Mul(- muNF, Add(SSD_W22_00, DSD_phiphi_00), W133NF.dummy),
                               Mul(2, sqrt(muNF), SSD_phi_00, Add(W02NF.dummy, W22NF.dummy)),
                               Mul(2, muNF, SSD_W133_00, W22NF.dummy),
                               Mul(muNF, Add(SSD_W034_00, Mul(2, DSD_phiW133_00)),
                                   phiNF.dummy),
                               Mul(2, sqrt(muNF), Add(SSD_W02_00, DSD_phiphi_00),
                                   phiNF.dummy)).subs(eq_subs).subs(extraparvals)

alpha4 = psiNF.dummy.dot(negativeRHS.actualcoord)

SS_W13_10 = first_order(1, 0, W13NF)
SS_W13_01 = first_order(0, 1, W13NF)
SS_phi_20 = first_order(2, 0, phiNF)
SS_phi_11 = first_order(1, 1, phiNF)
SS_phi_02 = first_order(0, 2, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[2], SS_W13_10), Mul(bNF[2], SS_W13_01),
                              Mul(aNF[4], SS_phi_10), Mul(bNF[4], SS_phi_01),                       
                              Mul(Pow(aNF[2], 2), SS_phi_20),
                              Mul(aNF[2], bNF[2], SS_phi_11), Mul(Pow(bNF[2], 2), SS_phi_02),
                              Mul(- muNF, Add(Mul(aNF[2], diffmatrix10),
                                              Mul(bNF[2], diffmatrix01)), W13NF.dummy),
                              Mul(- muNF, Add(Mul(aNF[4], diffmatrix10),
                                              Mul(bNF[4], diffmatrix01),
                                              Mul(Pow(aNF[2], 2), diffmatrix20),
                                              Mul(aNF[2], bNF[2], diffmatrix11),
                                              Mul(Pow(bNF[2], 2), diffmatrix02)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

alpha5 = psiNF.dummy.dot(negativeRHS.actualcoord)

SS_W123_10 = first_order(1, 0, W123NF)
SS_W123_01 = first_order(0, 1, W123NF)
DS_phiW04_00 = second_order(0, 0, phiNF, W04NF)
DS_phiW24_00 = second_order(0, 0, phiNF, W24NF)
DS_W13W02_00 = second_order(0, 0, W13NF, W02NF)
DS_W13W22_00 = second_order(0, 0, W13NF, W22NF)
DS_phiW02_10 = second_order(1, 0, phiNF, W02NF)
DS_phiW22_10 = second_order(1, 0, phiNF, W22NF)
DS_phiW02_01 = second_order(0, 1, phiNF, W02NF)
DS_phiW22_01 = second_order(0, 1, phiNF, W22NF)
TS_phiphiW13_00 = third_order(0, 0, phiNF, phiNF, W13NF)
TS_phiphiphi_10 = third_order(1, 0, phiNF, phiNF, phiNF)
TS_phiphiphi_01 = third_order(0, 1, phiNF, phiNF, phiNF)

SSD_W04_00 = first_order_diff(0, 0, W04NF)
SSD_W24_00 = first_order_diff(0, 0, W24NF)
SSD_W02_10 = first_order_diff(1, 0, W02NF)
SSD_W22_10 = first_order_diff(1, 0, W22NF)
SSD_W02_01 = first_order_diff(0, 1, W02NF)
SSD_W22_01 = first_order_diff(0, 1, W22NF)
DSD_phiW13_00 = second_order_diff(0, 0, phiNF, W13NF)
DSD_phiphi_10 = second_order_diff(1, 0, phiNF, phiNF)
DSD_phiphi_01 = second_order_diff(0, 1, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[2], SS_W123_10), Mul(bNF[2], SS_W123_01),
                              Mul(2, Add(Mul(2, DS_phiW04_00), DS_phiW24_00)),
                              Mul(2, Add(Mul(2, DS_W13W02_00), DS_W13W22_00)),
                              Mul(2, aNF[2], Add(Mul(2, DS_phiW02_10), DS_phiW22_10)),
                              Mul(2, bNF[2], Add(Mul(2, DS_phiW02_01), DS_phiW22_01)),
                              Mul(9, TS_phiphiW13_00),
                              Mul(3, aNF[2], TS_phiphiphi_10),
                              Mul(3, bNF[2], TS_phiphiphi_01),
                              Mul(- 2, muNF, SSD_phi_00, W24NF.dummy),
                              Mul(- muNF, Add(Mul(2, SSD_W02_00), - SSD_W22_00,
                                              DSD_phiphi_00), W13NF.dummy),
                              Mul(- muNF, Add(Mul(aNF[2], diffmatrix10),
                                              Mul(bNF[2], diffmatrix01)),
                                  W123NF.dummy),
                              Mul(- 2, muNF, Add(SSD_W13_00,
                                                 Mul(aNF[2], SSD_phi_10),
                                                 Mul(bNF[2], SSD_phi_01)),
                                  W22NF.dummy),
                              Mul(- muNF, Add(Mul(2, SSD_W04_00), - SSD_W24_00,
                                            Mul(aNF[2], Add(Mul(2, SSD_W02_10), - SSD_W22_10)),
                                            Mul(bNF[2], Add(Mul(2, SSD_W02_01), - SSD_W22_01)),
                                            Mul(2, DSD_phiW13_00),
                                            Mul(aNF[2], DSD_phiphi_10),
                                            Mul(bNF[2], DSD_phiphi_01)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

alpha6 = psiNF.dummy.dot(negativeRHS.actualcoord)

DS_phiW024_00 = second_order(0, 0, phiNF, W024NF)
DS_phiW224_00 = second_order(0, 0, phiNF, W224NF)
DS_W22W33_00 = second_order(0, 0, W22NF, W33NF)
DS_W123W02_00 = second_order(0, 0, W123NF, W02NF)
DS_W123W22_00 = second_order(0, 0, W123NF, W22NF)
TS_phiphiW123_00 = third_order(0, 0, phiNF, phiNF, W123NF)
TS_phiphiW33_00 = third_order(0, 0, phiNF, phiNF, W33NF)
TS_phiW02W02_00 = third_order(0, 0, phiNF, W02NF, W02NF)
TS_phiW22W02_00 = third_order(0, 0, phiNF, W22NF, W02NF)
TS_phiW22W22_00 = third_order(0, 0, phiNF, W22NF, W22NF)
Q4S_phiphiphiW02_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W22NF)
Q5S_phiphiphiphiphi_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, phiNF)

SSD_W024_00 = first_order_diff(0, 0, W024NF)
SSD_W224_00 = first_order_diff(0, 0, W224NF)
DSD_phiW123_00 = second_order_diff(0, 0, phiNF, W123NF)
DSD_phiW33_00 = second_order_diff(0, 0, phiNF, W33NF)
DSD_W02W02_00 = second_order_diff(0, 0, W02NF, W02NF)
DSD_W02W22_00 = second_order_diff(0, 0, W02NF, W22NF)
DSD_W22W22_00 = second_order_diff(0, 0, W22NF, W22NF)
TSD_phiphiW02_00 = third_order_diff(0, 0, phiNF, phiNF, W02NF)
Q4SD_phiphiphiphi_00 = fourth_order_diff(0, 0, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(2, Add(Mul(2, DS_phiW024_00), DS_phiW224_00)),
                              Mul(2, DS_W22W33_00), Mul(2, Add(Mul(2, DS_W123W02_00), DS_W123W22_00)),
                              Mul(3, Add(Mul(3, TS_phiphiW123_00), TS_phiphiW33_00)),
                              Mul(12, TS_phiW02W02_00),
                              Mul(6, Add(Mul(2, TS_phiW22W02_00), TS_phiW22W22_00)),
                              Mul(8, Add(Mul(3, Q4S_phiphiphiW02_00), Mul(2, Q4S_phiphiphiW22_00))),
                              Mul(10, Q5S_phiphiphiphiphi_00),
                              Mul(- 2, muNF, SSD_phi_00, W224NF.dummy),
                              Mul(- muNF, Add(Mul(2, SSD_W02_00),
                                              - SSD_W22_00, DSD_phiphi_00),
                                  W123NF.dummy),
                              Mul(- 3, muNF, Add(SSD_W22_00, DSD_phiphi_00),
                                  W33NF.dummy),
                              Mul(- 2, muNF, Add(SSD_W123_00, - SSD_W33_00,
                                                 Mul(4, DSD_phiW02_00),
                                                 Mul(2, TSD_phiphiphi_00)),
                                  W22NF.dummy),
                              Mul(- muNF, Add(Mul(2, SSD_W024_00),
                                              - SSD_W224_00,
                                              Mul(2, Add(DSD_phiW123_00, - DSD_phiW33_00)),
                                              Mul(4, Add(DSD_W02W02_00, - DSD_W02W22_00)),
                                              Mul(2, DSD_W22W22_00),
                                              Mul(6, TSD_phiphiW02_00),
                                              Mul(2, Q4SD_phiphiphiphi_00)),
                                  phiNF.dummy)).subs(eq_subs).subs(extraparvals)

alpha7 = psiNF.dummy.dot(negativeRHS.actualcoord)

if not os.path.exists('alphas'):
    os.makedirs('alphas')

os.chdir(os.getcwd() + '\\alphas')

for counter in range(1, 8):
    if counter<3:
        file = open(f'upsilon{counter}.txt', 'w')
        exec(f'file.write(mathematica_code(upsilon{counter}))')
        file.close()
    file = open(f'alpha{counter}.txt', 'w')
    exec(f'file.write(mathematica_code(alpha{counter}))')
    file.close()
    
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.getcwd() + '\\' + modelname)

file = open('Vectors.txt', 'w')
write_vector(phiNF, file)
write_vector(psiNF, file)
write_vector(W02NF, file)
write_vector(W22NF, file)
write_vector(W13NF, file)
write_vector(W123NF, file)
write_vector(W133NF, file)
write_vector(W33NF, file)
write_vector(W04NF, file)
write_vector(W024NF, file)
write_vector(W034NF, file)
write_vector(W24NF, file)
write_vector(W224NF, file)
write_vector(W234NF, file)
file.close()

try:
    eval(parameters_on_axes[0])
    eval(parameters_on_axes[0])
    file = open('Parameters on axes.txt', 'w')
    file.write(parameters_on_axes[0] + ',' + mathematica_code(intervalx[0]) + ',' + mathematica_code(intervalx[1]) + '\n')
    file.write(parameters_on_axes[1] + ',' + mathematica_code(intervaly[0]) + ',' + mathematica_code(intervaly[1]) + '\n')
    file.close()
except:
    print('The parameters on axes provided are not parameters of the system.')

file = open('Initial conditions for Turing bifurcation curves.txt','w')

for key in lines_to_search.keys():
    if isinstance(lines_to_search[key],list):
        for initialsolnum in range(len(lines_to_search[key])):
            if isfloat(str(lines_to_search[key][initialsolnum])):
                file.write(mathematica_code(eval(key)) + ',' + 
                           mathematica_code(eval(str(lines_to_search[key][initialsolnum]))) + '\n')
    else:
        if isfloat(str(lines_to_search[key])):
            file.write(mathematica_code(eval(key)) + ',' + mathematica_code(eval(str(lines_to_search[key]))) + '\n')
        
file.close()

file = open('Fixed parameter values.txt', 'w')
for key in parameters.keys():
    if key not in parameters_on_axes:
        if key not in parameter_functions.keys():
            file.write(mathematica_code(key) + ',' + mathematica_code(parameters[key]) + '\n')
        else:
            file.write(mathematica_code(key) + ',' + mathematica_code(parameter_functions[key]) + '\n')

file.close()

try:
    if len(names_of_parameters)==0:
        for parnum in range(2):
            names_of_parameters[parnum] = parameters_on_axes[parnum]
except:
    names_of_parameters = parameters_on_axes
    for parnum in range(2):
        names_of_parameters[parnum] = names_of_parameters[parnum]
    
file = open('Actual names of parameters.txt', 'w')
file.write(names_of_parameters[0] + ',' + names_of_parameters[1])
file.close()

file = open('Expansion parameters.txt', 'w')
file.write(mathematica_code(par1) + ',')
file.write(mathematica_code(par2))
file.close()

if not os.path.isfile('Plotter.nb'):
    shutil.copyfile(os.path.dirname(os.path.realpath(__file__)) + '\\Plotter.nb', 'Plotter.nb')
    
if modelname=='Nonlinear_SH':
    mu = simplify(Mul(Add(Pow(alpha2, 2), Mul(4, alpha1, alpha5)),
             Pow(Mul(4, Pow(alpha1, 2)), - 1)).subs(W13NF_eval).subs(W133NF_eval).\
                  subs(phiNF_eval).subs(psiNF_eval).subs(extraparvals))
    nalpha3 = simplify(Mul(Add(Mul(alpha2, Add(alpha3, - alpha4)),
                      Mul(2, alpha1, alpha6)),
                  Pow(Mul(2, Pow(alpha1, 2)), - 1)).subs(W04NF_eval).subs(W24NF_eval).\
                       subs(W13NF_eval).subs(W123NF_eval).\
                           subs(W133NF_eval).subs(phiNF_eval).\
                               subs(psiNF_eval))
    nalpha5 = simplify(Mul(alpha7, Pow(alpha1, - 1)).subs(W024NF_eval).subs(W224NF_eval).\
                       subs(W123NF_eval).subs(W133NF_eval).subs(W33NF_eval).\
                       subs(W02NF_eval).subs(W22NF_eval).subs(phiNF_eval).\
                           subs(psiNF_eval))
    nbeta1 = simplify(Mul(alpha3, Pow(alpha1, - 1)).subs(W034NF_eval).subs(W234NF_eval).\
                      subs(W123NF_eval).subs(W133NF_eval).subs(W02NF_eval).\
                          subs(W22NF_eval).subs(phiNF_eval).subs(psiNF_eval))
    nbeta2 = simplify(Mul(alpha4, Pow(alpha1, - 1)).subs(W034NF_eval).subs(W123NF_eval).\
                      subs(W133NF_eval).subs(W02NF_eval).subs(W22NF_eval).\
                      subs(phiNF_eval).subs(psiNF_eval))
    
    BNF = Add(nbeta2, - nbeta1)
    HNF = Mul(Add(nbeta2, nbeta1), Pow(4, - 1))
    
    hNF = simplify(Add(Pow(HNF, 2), Mul(BNF, HNF), - nalpha5))
    
    stab = simplify(Add(Mul(2, nalpha5), Mul(HNF, Add(BNF, Mul(4, HNF)))))