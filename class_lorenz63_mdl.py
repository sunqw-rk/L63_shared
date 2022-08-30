# Class S-I1-I2-H-R model
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pickle
import math
import pylab as pl
import matplotlib.ticker as mticker
sigma = 10.0
beta = 8.0/3.0
rho = 28.0

def lorenz63(t,x): 
    d = np.zeros(3)
    d[0] = sigma * (x[1] - x[0])
    d[1] = x[0] * (rho - x[2]) - x[1]
    d[2] = x[0] * x[1] - beta * x[2]

    # Return the state derivatives
    return d
    

class LORENZ63:
    def __init__(self, sigma, beta, rho):
        self.sigma = sigma
        self.beta = beta
        self.rho = rho
        self.params = [self.sigma, self.beta, self.rho]

    #--------------------------------------
    #using odeint to slove the siir system
    #--------------------------------------
    def lorenz63(self,x): 
        d = np.zeros(3)
        d[0] = sigma * (x[1] - x[0])
        d[1] = x[0] * (rho - x[2]) - x[1]
        d[2] = x[0] * x[1] - beta * x[2]

    # Return the state derivatives
        return d
    
    def run_rk4(self, time_start, time_end, dt, state_all0):
        n = np.arange(time_start, time_end+1e-5, dt)
        y = state_all0
        y = np.array(y)
        stacked_array = y
        for i in range(len(n)-1):
            dx1 = np.array(self.lorenz63(y))
            x1 = y+ dx1 * dt * 0.5
            dx2 = np.array(self.lorenz63(x1))
            x2 = y + dx2 * dt * 0.5
            dx3 = np.array(self.lorenz63(x2))
            x3 = y + dx3 * dt
            dx4 = np.array(self.lorenz63(x3))
            y = y + (dt / 6.0) *(dx1 + 2.0 * dx2 + 2.0 * dx3 + dx4)
            time_start = time_start + dt
            stacked_array = np.vstack((stacked_array, y))
        return stacked_array
    
    def run(self, x0, t, t_output):
        x = solve_ivp(lorenz63, (t[0],t[-1]), x0, t_eval= t_output) 
        return x