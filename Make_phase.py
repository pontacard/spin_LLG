import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin.one_spin_LLG import A_spin
from one_spin.STO_gif import STO_gif

S0 = [0, 0, 1]

t = [0, 500]  # t(時間)が0〜100まで動き、その時のfを求める。
t_eval = np.linspace(*t, 70000)

plotB = [[0, 0, -1.2], [0, 0, 2.4]]

gamma = 0.17
mu_0 = 1.26
mu_h_div_2e = [0.824, -21]
sta_M = [1.2, 0]  # 飽和磁化(T)で入れる
jac = [33.44 * 4.02, 10]
d = [2, -9]
Hacn = mu_h_div_2e[0] * jac[0] / (sta_M[0] * d[0])
Haco = mu_h_div_2e[1] + jac[1] - (sta_M[1] + d[1])
Hac = Hacn * (10 ** Haco) * 1000 / gamma  # 最後の1000はmTにするため
print(Hac)

jdc = [2.2 * 4, 10]
Hdcn = mu_h_div_2e[0] * jdc[0] / (sta_M[0] * d[0])
Hdco = mu_h_div_2e[1] + jdc[1] - (sta_M[1] + d[1])
Hdc = Hacn * (10 ** Haco) * 1000 / gamma  # 最後の1000はmTにするため

Kx = 0
Ky = 0
Kz = 1481 - 1448

for jac in [120,130,133,134,134.43,134.8,135,140]:
    for omega in [1.5,2,2.5,3,3.5,4]:
        spin = STO_gif(0.005, gamma, [0, 0, mu_0 * 159], S0, t, t_eval, [0, 0, 0], [2, 2, 0], [0, 0, 0], Kx, Ky,
                    mu_0 * Kz, 0, 0, 40, [0, Hac, 0], [0, Hdc, 0], 0.288, 0.537)
spin.make_gif()