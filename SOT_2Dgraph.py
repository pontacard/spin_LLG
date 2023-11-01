import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from SOT import SOT
class SOT_2D(SOT):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,spin_flow,Kx,Ky,Kz,beta,start,stop):
        super().__init__(alpha,gamma,B,S0,t,t_eval,spin_flow,Kx,Ky,Kz,beta,start,stop)


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
        #plt.show()
        plt.plot(t, y, label = "Sy")
        #plt.savefig("reverse_y.png")
        #plt.show()
        plt.plot(t, z,  label = "Sz")
        plt.ylim(-1.1,1.1)

        plt.xlabel("time(ns)")
        plt.ylabel("S")
        #plt.axhline(color = 'k')
        plt.legend()

        plt.savefig("TypeX_0-4_1200step_xyzS.pdf")
        plt.show()


if __name__ == '__main__':
    S0 = [1, 0, 0]

    t = [0, 4]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 1200)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 2
    mu_h_div_2e = [0.824, -21]
    sta_M = [1.4, 0]  # 飽和磁化(T)で入れる
    theta = [-2, -1]
    j = [2.5, 11]
    d = [1, -9]
    Hsn = mu_h_div_2e[0] * theta[0] * j[0] / (sta_M[0] * d[0])
    Hso = mu_h_div_2e[1] + theta[1] + j[1] - (sta_M[1] + d[1])
    Hs = Hsn * (10 ** Hso) * 1000 / gamma  # 最後の1000はmTにするため
    print(Hs)

    spin = SOT_2D(0.02, gamma, [0, 0, 7], S0, t, t_eval, [0, -Hs, 0], 2, 0, -50, 0, 0, 1)
    spin.get_graph()