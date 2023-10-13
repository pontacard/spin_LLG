import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from SOT import SOT
class SOT_2D(SOT):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,spin_flow,K,sc,fi,start,stop):
        super().__init__(alpha,gamma,B,S0,t,t_eval,spin_flow,K,sc,fi,start,stop)


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
        plt.plot(t, S, label = "|S|", ls = "--")
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

        plt.savefig("TypeY_0-0.03_4000step_xyzS.pdf")
        plt.show()


if __name__ == '__main__':
    S0 = [0, 0, -1]

    t = [0, 1.2]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 4000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    spin = SOT_2D(0.1, 28, [0.01, 0, 0], S0, t, t_eval, [0,0, -60], 84,0.1,1,0,0.03)
    spin.get_graph()