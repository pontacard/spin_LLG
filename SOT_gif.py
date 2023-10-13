import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin_gif import one_spin_gif

class SOT(one_spin_gif):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,spin_flow,K,sc,fi,start,stop):
        super().__init__(alpha,gamma,B,S0,t,t_eval,plotB)
        self.spin_flow = spin_flow
        self.K = K
        self.sc = sc
        self.fi = fi
        self.start = start
        self.stop = stop

    def func_S(self,t,S):  # 関数を設定
        if t>= self.start and t <= self.stop:
            B_sc = [self.B[0] + self.sc * self.spin_flow[0], self.B[1] + self.sc * self.spin_flow[1],
                    self.B[2] + self.sc * self.spin_flow[2]]
            B_fi = [self.B[0] + self.fi * self.spin_flow[0], self.B[1] + self.fi * self.spin_flow[1],
                    self.B[2] + self.fi * self.spin_flow[2]]
        else:
            B_sc = [self.B[0], self.B[1], self.B[2]]
            B_fi = [self.B[0], self.B[1], self.B[2]]

        Snorm = np.linalg.norm(S)
        KneZ = [self.K * S[1] * S[2] / (Snorm * Snorm), -self.K * S[2] * S[0] / (Snorm * Snorm), 0]

        dSxdt = - self.gamma * (B_sc[1] * S[2] - B_sc[2] * S[1]) + KneZ[0] - (self.alpha / Snorm) * (
                    S[1] * (self.gamma * (B_fi[0] * S[1] - B_fi[1] * S[0]) + KneZ[2]) - S[2] * (
                        self.gamma * (B_fi[2] * S[0] - B_fi[0] * S[2]) + KneZ[1]))
        dSydt = - self.gamma * (B_sc[2] * S[0] - B_sc[0] * S[2]) + KneZ[1] - (self.alpha / Snorm) * (
                    S[2] * (self.gamma * (B_fi[1] * S[2] - B_fi[2] * S[1]) + KneZ[0]) - S[0] * (
                        self.gamma * (B_fi[0] * S[1] - B_fi[1] * S[0]) + KneZ[2]))
        dSzdt = - self.gamma * (B_sc[0] * S[1] - B_sc[1] * S[0]) + KneZ[2] - (self.alpha / Snorm) * (
                    S[0] * (self.gamma * (B_fi[2] * S[0] - B_fi[0] * S[2]) + KneZ[1]) - S[1] * (
                        self.gamma * (B_fi[1] * S[2] - B_fi[2] * S[1]) + KneZ[0]))
        dSdt = [dSxdt, dSydt, dSzdt]
        # print(dSdt)

        return dSdt

    def make_gif(self):
        self.fig, self.ax = plt.subplots(subplot_kw=dict(projection="3d"))

        self.Sol = sc.integrate.solve_ivp(self.func_S, self.t, self.S0, t_eval=self.t_eval)
        self.S = self.Sol.y

        self.quiveraa = self.ax.quiver(*self.get_spin_vec(0))

        r  = np.linalg.norm(self.S0) # 半径を指定
        theta_1_0 = np.linspace(0, np.pi * 2, 100)  # θ_1は&#91;0,π/2]の値をとる
        theta_2_0 = np.linspace(0, np.pi, 100)  # θ_2は&#91;0,π/2]の値をとる
        theta_1, theta_2 = np.meshgrid(theta_1_0, theta_2_0)  # ２次元配列に変換
        x = np.cos(theta_2) * np.sin(theta_1) * r  # xの極座標表示
        y = np.sin(theta_2) * np.sin(theta_1) * r  # yの極座標表示
        z = np.cos(theta_1) * r  # zの極座標表示

        self.ax.plot_surface(x, y, z, alpha=0.2)  # 球を３次元空間に表示

        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_zlabel("z")
        self.ax.text(0.12, 0.15, -0.3, '-B',
                size=20, )

        self.ax.view_init(elev=20, azim=20)
        self.ax.set_xlim(-1.2, 1.2)
        self.ax.set_ylim(-1.2, 1.2)
        self.ax.set_zlim(-1.2, 1.2)
        self.ax.set_aspect('equal')

        ani = FuncAnimation(self.fig, self.update, frames=len(self.Sol.t), interval=1)
        #ani.save("SOT_Type_Z.gif",writer='imagemagick')
        plt.show()


if __name__ == '__main__':
    S0 = [0, 0, -1]

    t = [0, 1]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 4000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    spin = SOT(0.1, 28, [0.01, 0, 0], S0, t, t_eval, [0, -60, 0], 84,0.1,1,0,0.01)
    spin.make_gif()
