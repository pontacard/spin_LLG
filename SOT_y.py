import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from SOT_gif import SOT
from SOT_2Dgraph import SOT_2D


S0 = [np.sqrt(0.000001), -np.sqrt(0.999999), 0]

t = [0,1.5] # t(時間)が0〜100まで動き、その時のfを求める。
t_eval = np.linspace(*t, 500)

plotB = [[0, 0, -1.2], [0, 0, 2.4]]

spin = SOT(0.1,-28,[0,-4,0],S0,t,t_eval, plotB,[0,-40,0])
spin.make_gif()

D2_spim = SOT_2D(0.1,-28,[0,-4,0],S0,t,t_eval,[0,-40,0])
D2_spim.get_graph()

