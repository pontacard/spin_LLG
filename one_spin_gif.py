import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin_LLG import A_spin

class one_spin_gif(A_spin):
    def __init__(self,alpha,gamma,B,S0,t,plotB):
        super().__init__(alpha,gamma,B,S0,t)
        self.plotB = plotB

    def get_spin_vec(self,t):
        Ox, Oy, Oz = 0, 0, 0
        x = self.S[t][0]
        y = self.S[t][1]
        z = self.S[t][2]
        self.ax.plot(x, y, z, marker='o', markersize=2.5, color='b')
        return Ox, Oy, Oz, x, y, z

    def make_gif(self):
        self.fig, self.ax = plt.subplots(subplot_kw=dict(projection="3d"))

        self.S = sc.integrate.odeint(self.func_S, self.S0, self.t)

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

        ani = FuncAnimation(self.fig, self.update, frames=len(self.t), interval=1)
        # ani.save("reverse_spin.gif",writer='imagemagick')
        plt.show()

if __name__ == '__main__':
    S0 = [0.1, 0, -0.9]

    t = np.linspace(0, 150, 15000)  # t(時間)が0〜100まで動き、その時のfを求める。

    plotB = [[0,0,-1.2],[0,0,2.4]]

    spin = one_spin_gif(0, 1, [0, 0, 5], S0, t,plotB)
    spin.make_gif()