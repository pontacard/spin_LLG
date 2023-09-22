import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin_gif import one_spin_gif

class SOT(one_spin_gif):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,plotB,spin_flow):
        super().__init__(alpha,gamma,B,S0,t,t_eval,plotB)
        self.spin_flow = spin_flow

    def func_S(self,t,S):  # 関数を設定
        B = self.B + self.spin_flow
        Snorm = np.linalg.norm(S)
        dSxdt = - self.gamma * (B[1] * S[2] - B[2] * S[1]) - (self.gamma * self.alpha / Snorm) * (
                    S[1] * (S[0] * (B[1] + self.spin_flow[1]) - S[1] * B[0]) - S[2] * (S[2] * B[0] - S[0] * B[2]))
        dSydt = - self.gamma * (B[2] * S[0] - B[0] * S[2]) - (self.gamma * self.alpha / Snorm) * (
                    S[2] * (S[1] * B[2] - S[2] * B[1]) - S[0] * (S[0] * (B[1] + self.spin_flow[1]) - S[1] * B[0]))
        dSzdt = - self.gamma * (B[0] * S[1] - B[1] * S[0]) - (self.gamma * self.alpha / Snorm) * (
                    S[0] * (S[2] * B[0] - S[0] * B[2]) - S[1] * (S[1] * B[2] - S[2] * (B[1] + self.spin_flow[1])))


        dSdt = [dSxdt, dSydt, dSzdt]

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

        self.ax.quiver(self.plotB[0][0], self.plotB[0][1], self.plotB[0][2], self.plotB[1][0], self.plotB[1][1], self.plotB[1][2], color='k', arrow_length_ratio=0.18, linewidth=3, label='B')
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
        #ani.save("SOT_Type_Y.gif",writer='imagemagick')
        plt.show()


if __name__ == '__main__':
    S0 = [0, 0, -1]

    t = [0,1] # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 1000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    spin = SOT(0.1,-28,[0.01,0,-4],S0,t,t_eval, plotB,[0,-10,0])
    spin.make_gif()