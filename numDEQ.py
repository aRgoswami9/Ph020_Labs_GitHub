from itertools import islice
import matplotlib.pyplot as plt
import numpy as np
import math

def xv_gen_explicit(x_0, v_0, h):
        '''Generates x,v-values using the explicit Euler method'''
        xv_tuple = (x_0, v_0)
        while True:
            x, v = xv_tuple
            yield xv_tuple
            xv_tuple = (x + h * v, v - h * x)

def xv_gen_implicit(x_0, v_0, h):
        '''Generates x,v-values using the implicit Euler method'''
        xv_tuple = (x_0, v_0)
        while True:
            x, v = xv_tuple
            yield xv_tuple
            xv_tuple = ((x + h * v) / (1 + h**2), \
                (v - h * x) / (1 + h**2))

def xv_gen_symplectic(x_0, v_0, h):
        '''Generates x,v-values using the symplectic Euler method'''
        xv_tuple = (x_0, v_0)
        while True:
            x, v = xv_tuple
            yield xv_tuple
            xv_tuple = (x + h * v, -h * x + v - h * h * v)

def springEuler(x_0, v_0, h, numSteps):
    '''Numerically investigates the motion of a mass on 
    a string, using the explicit Euler method. Takes initial
    conditions x_0, v_0, step size h, and total number of steps 
    as parameters. Plots X and V versus time'''
    xArr = np.zeros(numSteps)
    vArr = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_explicit(x_0, v_0, h), numSteps)):
        xArr[i], vArr[i] = e
    tArr = h * np.arange(numSteps)
    fig = plt.figure()
    plt.plot(tArr, xArr, 'k-', tArr, vArr, 'b-')
    plt.savefig('explicit_spring_motion.pdf')

def plotError(x_0, v_0, h, numSteps): 
    '''Given the initial position, velocity, step size, and 
    number of steps as parameters, plots the error as a function of
    time between the numerical solution and the analytical solution.
    Returns the maximum error between the numerical and analytical
    solution'''
    xArr = np.zeros(numSteps)
    vArr = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_explicit(x_0, v_0, h), numSteps)):
        xArr[i], vArr[i] = e
    tArr = h * np.arange(numSteps)
    xArrAnal = x_0 * np.cos(tArr) + v_0 * np.sin(tArr)
    vArrAnal = -x_0 * np.sin(tArr) + v_0 * np.cos(tArr)
    fig1 = plt.figure()
    plt.plot(tArr, np.absolute(xArrAnal - xArr), 'k-', \
             tArr, np.absolute(vArrAnal - vArr), 'b-')
    plt.savefig('explicit_spring_error.pdf')
    return np.max(np.absolute(xArrAnal - xArr))

def plotTruncationError(x_0, v_0, h, finalTime):
    '''Shows that the truncation error is roughly proportional
    to h for small values of h'''
    hList = [h, h / 2, h / 4, h / 8, h / 16]
    ErrList = []
    for i in range(len(hList)):
        numSteps = int(finalTime / hList[i])
        ErrList.append(plotError(x_0, v_0, hList[i], numSteps))
    fig = plt.figure()
    plt.plot(hList, ErrList, 'ko-')
    plt.savefig('truncation_error.pdf')

def plotEnergyEvolution(x_0, v_0, h, numSteps):
    '''Plots the energy evolution of both the numerical spring 
    system (in black) and the analytical spring system (in blue)
    with given parameters'''
    xArr = np.zeros(numSteps)
    vArr = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_explicit(x_0, v_0, h), numSteps)):
        xArr[i], vArr[i] = e
    tArr = h * np.arange(numSteps)
    xArrAnal = x_0 * np.cos(tArr) + v_0 * np.sin(tArr)
    vArrAnal = -x_0 * np.sin(tArr) + v_0 * np.cos(tArr)
    fig = plt.figure()
    plt.plot(tArr, xArr**2 + vArr**2, 'k-', \
        tArr, xArrAnal**2 + vArrAnal**2, 'b-')
    plt.savefig('energy_evolution_explicit.pdf')

def implicitEuler(x_0, v_0, h, numSteps):
    '''Investigates the spring system with the implicit Euler method 
    (black),compares energy evolution and global errors with explicit 
    Euler method (blue)'''
    tArr = h * np.arange(numSteps)
    xArrImp = np.zeros(numSteps)
    vArrImp = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_implicit(x_0, v_0, h), numSteps)):
        xArrImp[i], vArrImp[i] = e
    xArrExp = np.zeros(numSteps)
    vArrExp = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_explicit(x_0, v_0, h), numSteps)):
        xArrExp[i], vArrExp[i] = e
    xArrAnal = x_0 * np.cos(tArr) + v_0 * np.sin(tArr)
    vArrAnal = -x_0 * np.sin(tArr) + v_0 * np.cos(tArr)
    figEnergy = plt.figure()
    plt.plot(tArr, xArrImp**2+vArrImp**2, 'k-', \
        tArr, xArrExp**2+vArrExp**2, 'b-')
    plt.savefig('energy_evolution_implicit_explicit.pdf')
    figError = plt.figure()
    plt.plot(tArr, np.absolute(xArrImp - xArrAnal), 'k-', \
        tArr, np.absolute(xArrExp - xArrAnal), 'b-')
    plt.savefig('implicit_global_error.pdf')

def phaseSpace(x_0, v_0, h, numSteps):
    '''Investigates the phase space geometry of the 
    implicit (black),  explicit (blue), and symplectic
    (red) Euler methods'''
    xArrImp = np.zeros(numSteps)
    vArrImp = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_implicit(x_0, v_0, h), numSteps)):
        xArrImp[i], vArrImp[i] = e
    xArrExp = np.zeros(numSteps)
    vArrExp = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_explicit(x_0, v_0, h), numSteps)):
        xArrExp[i], vArrExp[i] = e
    xArrSymp = np.zeros(numSteps)
    vArrSymp = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_symplectic(x_0, v_0, h), numSteps)):
        xArrSymp[i], vArrSymp[i] = e
    fig = plt.figure()
    plt.plot(xArrImp, vArrImp, 'k-')
    plt.savefig('implicit_phasespace.pdf')
    fig2 = plt.figure()
    plt.plot(xArrExp, vArrExp, 'b-')
    plt.savefig('explicit_phasespace.pdf')
    fig3 = plt.figure()
    plt.plot(xArrSymp, vArrSymp, 'r-')
    plt.savefig('symplectic_phasespace.pdf')

def plotSymplecticEnergyEvolution(x_0, v_0, h, numSteps):
    '''Plots the energy evolution of symplectic Euler method
    for spring with given parameters'''
    xArrSymp = np.zeros(numSteps)
    vArrSymp = np.zeros(numSteps)
    for i, e in enumerate(islice(xv_gen_symplectic(x_0, v_0, h), numSteps)):
        xArrSymp[i], vArrSymp[i] = e
    tArr = h * np.arange(numSteps)
    xArrAnal = x_0 * np.cos(tArr) + v_0 * np.sin(tArr)
    vArrAnal = -x_0 * np.sin(tArr) + v_0 * np.cos(tArr)
    fig = plt.figure()
    plt.plot(tArr, xArrSymp**2 + vArrSymp**2, 'k-', \
        tArr, xArrAnal**2 + vArrAnal**2, 'b-')
    plt.savefig('symplectic_energy_evolution.pdf')

# comment added