import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin_LLG import A_spin

class SOT(A_spin):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,spin_flow,K,sc,fi,start,stop): #scはSlonczewski-like torque、fiはfield-like torqueをそれぞれ実行的磁場と考えた時のバイアス
        super().__init__(alpha,gamma,B,S0,t,t_eval)
        self.spin_flow = spin_flow
        self.K = K
        self.sc = sc
        self.fi =  fi
        self.start = start
        self.stop = stop


    def func_S(self,t,S):  # 関数を設定
        B_dm = []
        B_fi = []
        if t>= self.start and t <= self.stop:
            B_dm= [self.B[0] + (self.sc/self.alpha) * self.spin_flow[0] ,self.B[1] + (self.sc/self.alpha) *self.spin_flow[1], self.B[2] + (self.sc/self.alpha) * self.spin_flow[2]]
            B_fi= [self.B[0] + self.fi * self.spin_flow[0] ,self.B[1] + self.fi *self.spin_flow[1], self.B[2] + self.fi * self.spin_flow[2]]
        elif t>= 1.4 and t <= 1.5:
            B_dm = [self.B[0] + self.sc * self.spin_flow[0], self.B[1] + self.sc * self.spin_flow[1],
                    self.B[2] + self.sc * self.spin_flow[2]]
            B_fi = [self.B[0] + self.fi * self.spin_flow[0], self.B[1] + self.fi * self.spin_flow[1],
                    self.B[2] + self.fi * self.spin_flow[2]]

        else:
            B_dm = [self.B[0],self.B[1], self.B[2]]
            B_fi = [self.B[0],self.B[1], self.B[2]]

        Snorm = np.linalg.norm(S)
        KneZ = [self.K * S[1] * S[2]/(Snorm * Snorm), -self.K * S[2] * S[0]/(Snorm * Snorm), 0]

        dSxdt = - self.gamma * (B_fi[1] * S[2] - B_fi[2] * S[1]) + KneZ[0] - (self.alpha/Snorm) * (S[1] * (self.gamma * (B_dm[0] * S[1] - B_dm[1] * S[0]) + KneZ[2]) - S[2]* (self.gamma * (B_dm[2] * S[0] - B_dm[0] * S[2]) + KneZ[1]))
        dSydt = - self.gamma * (B_fi[2] * S[0] - B_fi[0] * S[2]) + KneZ[1] - (self.alpha/Snorm) * (S[2] * (self.gamma * (B_dm[1] * S[2] - B_dm[2] * S[1]) + KneZ[0]) - S[0]* (self.gamma * (B_dm[0] * S[1] - B_dm[1] * S[0]) + KneZ[2]))
        dSzdt = - self.gamma * (B_fi[0] * S[1] - B_fi[1] * S[0]) + KneZ[2] - (self.alpha/Snorm) * (S[0] * (self.gamma * (B_dm[2] * S[0] - B_dm[0] * S[2]) + KneZ[1]) - S[1]* (self.gamma * (B_dm[1] * S[2] - B_dm[2] * S[1]) + KneZ[0]))
        dSdt = [dSxdt, dSydt, dSzdt]
        #print(dSdt)

        return dSdt


if __name__ == '__main__':
    S0 = [0, 0, -1]

    t = [0, 0.5]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 300)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]
    sc = 1
    fi = 0

    spin = SOT(0.1, 28, [0.1, 0, 0], S0, t, t_eval, [0, 3, 0],40,sc,fi ,0,0.02)
    spin.doit()