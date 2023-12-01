import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin.one_spin_LLG import A_spin

class STO_spin(A_spin):
    def __init__(self,alpha,gamma,B,S0,t, t_eval,Amp,omega,theta,Kx,Ky,Kz,beta,start,stop,STO_effH,lambdaa,eta):
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
        self.STO_effH = STO_effH
        self.lambdaa = lambdaa
        self.eta = eta



    def func_S(self,t,S):  # 関数を設定
        Snorm = np.linalg.norm(S)
        Bk = np.array([self.Kx * S[0] / (Snorm * Snorm), self.Ky * S[1] / (Snorm * Snorm),self.Kz * S[2] / (Snorm * Snorm)])
        sot_torque = [0,0,0]
        p = np.array([0,1,0])
        SOT_angle = 0
        if t >= self.start and t <= self.stop:
            B = np.array(self.B) + Bk + np.array([self.Amp[0] * np.cos(self.omega[0] * t +  self.theta[0]) , self.Amp[1] * np.cos(self.omega[1] * t + self.theta[1]),0])
            #print(B)
            SOT_angle = self.eta/(1 + self.lambdaa * np.dot(S,p))
            sot_torque = [self.gamma  * (S[1] * (self.STO_effH[0] * S[1] - self.STO_effH[1] * S[0]) - S[2] * (
                    self.STO_effH[2] * S[0] - self.STO_effH[0] * S[2])), self.gamma * (
                                 S[2] * (self.STO_effH[1] * S[2] - self.STO_effH[2] * S[1]) - S[0] * (
                                 self.STO_effH[0] * S[1] - self.STO_effH[1] * S[0])), self.gamma * (
                                 S[0] * (self.STO_effH[2] * S[0] - self.STO_effH[0] * S[2]) - S[1] * (
                                 self.STO_effH[1] * S[2] - self.STO_effH[2] * S[1]))]
        else:
            B = np.array(self.B) + Bk

        Snorm = np.linalg.norm(S)
        dSxdt = - self.gamma * (B[1] * S[2] - B[2] * S[1]) - SOT_angle * sot_torque[0] - (self.gamma * self.alpha/Snorm) * (S[1] * (S[0] * B[1] - S[1] * B[0]) - S[2]* (S[2] * B[0] - S[0] * B[2]))
        dSydt = - self.gamma * (B[2] * S[0] - B[0] * S[2]) - SOT_angle * sot_torque[1] - (self.gamma * self.alpha/Snorm) * (S[2] * (S[1] * B[2] - S[2] * B[1]) - S[0]* (S[0] * B[1] - S[1] * B[0]))
        dSzdt = - self.gamma * (B[0] * S[1] - B[1] * S[0]) - SOT_angle * sot_torque[2] - (self.gamma * self.alpha/Snorm) * (S[0] * (S[2] * B[0] - S[0] * B[2]) - S[1]* (S[1] * B[2] - S[2] * B[1]))
        dSdt = [dSxdt, dSydt, dSzdt]

        return dSdt

if __name__ == '__main__':
    S0 = [1, 0, 0]

    t = [0, 4]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 1200)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 0.17
    mu_h_div_2e = [0.824, -21]
    sta_M = [1.4, 0]  # 飽和磁化(T)で入れる
    j = [2.5, 11]
    d = [2, -9]
    Hacn = mu_h_div_2e[0] *  j[0] / (sta_M[0] * d[0])
    Haco = mu_h_div_2e[1] + j[1] - (sta_M[1] + d[1])
    Hac = Hacn * (10 ** Haco) * 1000 / gamma  # 最後の1000はmTにするため
    print(Hac)

    spin = STO_spin(0.005, gamma, [0, 0, 200], S0, t, t_eval, [Hac, Hac, 0],[6.8,6.8,0],[0,0,0], 0, 0, 1862 - 12 * 1400, 0, 0, 0.5)
    spin.doit()