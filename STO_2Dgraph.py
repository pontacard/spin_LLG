import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from one_spin.STO import STO_spin
class one_spin_2D(STO_spin):
    def __init__(self, alpha, gamma, B, S0, t, t_eval, H_Amp, omega, theta, Kx, Ky, Kz, beta, start, stop, STO_ac_effH,
                 STO_dc_effH, lambdaa, eta):
        super().__init__(alpha, gamma, B, S0, t, t_eval, H_Amp, omega, theta, Kx, Ky, Kz, beta, start, stop,
                         STO_ac_effH, STO_dc_effH, lambdaa, eta)


    def get_graph(self):
        self.Sol = sc.integrate.solve_ivp(self.func_S, self.t, self.S0, t_eval=self.t_eval)
        self.S = self.Sol.y
        t = self.Sol.t
        #print(t)

        x = self.S[0]
        y = self.S[1]
        z = self.S[2]
        S = np.sqrt(x*x + y*y + z*z)
        #print(S)
        #plt.plot(t, S, label = "|S|", ls = "--")
        #plt.savefig("reverse_|S|.png")
        #plt.plot(t,x)
        plt.plot(t, x, label="Sx")
        #plt.savefig("reverse_x.png")
        plt.show()

        plt.plot(t, z, label="Sz")
        plt.show()

        N = len(x)  # サンプル数
        f_s = 400  # サンプリングレート f_s[Hz] (任意)
        dt = (self.t[1] - self.t[0]) / len(self.t_eval)  # サンプリング周期 dt[s]

        y_fft = np.fft.fft(x)  # 離散フーリエ変換
        freq = np.fft.fftfreq(N, d=dt)  # 周波数を割り当てる（※後述）
        Amp = abs(y_fft / (N / 2))  # 音の大きさ（振幅の大きさ）
        plt.plot(freq[1:int(N / 2)], Amp[1:int(N / 2)])  # A-f グラフのプロット
        # plt.xscale("log")  # 横軸を対数軸にセット
        plt.xlabel("GHz")

        plt.show()
        plt.plot(t, y, label = "Sy")
        #plt.savefig("reverse_y.png")
        #plt.show()
        #plt.plot(t, z,  label = "Sz")
        plt.ylim(-1.1,1.1)

        plt.xlabel("time(ns)")
        plt.ylabel("S")
        #plt.axhline(color = 'k')
        plt.legend()

        #plt.savefig("no_dump_xy.pdf")
        plt.show()



if __name__ == '__main__':
    S0 = [3/5, 0, 4/5]

    t = [0, 40]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 1000000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 0.17
    mu_0 = 1.26
    mu_h_div_2e = [0.824, -21]
    sta_M = [1.2, 0]  # 飽和磁化(T)で入れる
    jac = [3.5, 12]
    d = [2, -9]
    Hacn = mu_h_div_2e[0] * jac[0] / (sta_M[0] * d[0])
    Haco = mu_h_div_2e[1] + jac[1] - (sta_M[1] + d[1])
    Hac = Hacn * (10 ** Haco) * 1000 / gamma  # 最後の1000はmTにするため
    print(Hac)

    jdc = [2.2, 11]
    Hdcn = mu_h_div_2e[0] * jdc[0] / (sta_M[0] * d[0])
    Hdco = mu_h_div_2e[1] + jdc[1] - (sta_M[1] + d[1])
    Hdc = Hacn * (10 ** Haco) * 1000 / gamma  # 最後の1000はmTにするため

    Kx = 0
    Ky = 0
    Kz = 1481 - 1448

    spin = one_spin_2D(0.005, gamma, [0, 0, 200], S0, t, t_eval, [0, 0, 0],[5.2,5.2,0],[0,0,0], Kx,Ky, mu_0 * Kz, 0, 0, 4000,[0,Hac,0],[0,Hdc,0],0.288,0.537)
    print(spin.history())
    spin.get_graph()
