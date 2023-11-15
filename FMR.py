
import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin_LLG import A_spin

class FMR_spin(A_spin):
    def __init__(self,alpha,gamma,B,S0,t, t_eval,Amp,omega,theta,Kx,Ky,Kz,beta,start,stop):
        super().__init__(alpha, gamma, B, S0, t, t_eval)
        self.Amp = Amp
        self.omega = omega
        self.Kx = Kx
        self.Ky = Ky
        self.Kz = Kz
        self.beta = beta
        self.start = start
        self.stop = stop
        self.theta = theta

    def func_S(self,t,S):  # 関数を設定
        S = self.S
        Snorm = np.linalg.norm(S)
        Bk = np.array([self.Kx * S[0] / (Snorm * Snorm), self.Ky * S[1] / (Snorm * Snorm),self.Kz * S[2] / (Snorm * Snorm)])
        if t >= self.start and t <= self.stop:
            B = np.array(self.B) + Bk + np.array([self.Amp[0] * np.cos(self.omega[0] * t) , self.Amp[1] * np.sin(self.omega[1] * t),self.Amp[2] * np.sin(self.omega[2] * t)])
            #print(B)
        else:
            B = np.array(self.B) + Bk

        Snorm = np.linalg.norm(S)
        dSxdt = - self.gamma * (B[1] * S[2] - B[2] * S[1]) - (self.gamma * self.alpha/Snorm) * (S[1] * (S[0] * B[1] - S[1] * B[0]) - S[2]* (S[2] * B[0] - S[0] * B[2]))
        dSydt = - self.gamma * (B[2] * S[0] - B[0] * S[2]) - (self.gamma * self.alpha/Snorm) * (S[2] * (S[1] * B[2] - S[2] * B[1]) - S[0]* (S[0] * B[1] - S[1] * B[0]))
        dSzdt = - self.gamma * (B[0] * S[1] - B[1] * S[0]) - (self.gamma * self.alpha/Snorm) * (S[0] * (S[2] * B[0] - S[0] * B[2]) - S[1]* (S[1] * B[2] - S[2] * B[1]))
        dSdt = [dSxdt, dSydt, dSzdt]

        return dSdt

if __name__ == '__main__':
    S0 = [0, 0, 1]

    t = [0,15] # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 100)
    mu_0 = 1.2

    spin = FMR_spin(0.01, 0.17, [0,0,-5],S0,t,t_eval,[10,10,0],[0.8,0.8,0],[0,0,0],0 ,0 ,mu_0 * 100, 0, 0, 5,0)
    spin.doit()
