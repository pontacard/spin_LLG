import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin_LLG import A_spin

class SOT(A_spin):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,spin_flow):
        super().__init__(alpha,gamma,B,S0,t,t_eval)
        self.spin_flow = spin_flow

    def func_S(self,t,S):  # 関数を設定
        B = self.B

        Snorm = np.linalg.norm(S)
        dSxdt = - self.gamma * (B[1] * S[2] - B[2] * S[1]) - (self.gamma * self.alpha/Snorm) * (S[1] * (S[0] * (B[1] + self.spin_flow[1]) - S[1] * B[0]) - S[2]* (S[2] * B[0] - S[0] * B[2]))
        dSydt = - self.gamma * (B[2] * S[0] - B[0] * S[2]) - (self.gamma * self.alpha/Snorm) * (S[2] * (S[1] * B[2] - S[2] * B[1]) - S[0]* (S[0] * (B[1] + self.spin_flow[1]) - S[1] * B[0]))
        dSzdt = - self.gamma * (B[0] * S[1] - B[1] * S[0]) - (self.gamma * self.alpha/Snorm) * (S[0] * (S[2] * B[0] - S[0] * B[2]) - S[1]* (S[1] * B[2] - S[2] * (B[1] + self.spin_flow[1])))
        dSdt = [dSxdt, dSydt, dSzdt]
        print(dSdt)
        print(B)

        return dSdt


if __name__ == '__main__':
    S0 = [0, 0, -1]

    t = [0, 1]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 300)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    spin = SOT(0.1, -28, [0, 0, -4], S0, t, t_eval, [0, -10, 0])
    spin.doit()