import numpy as np

pulse = []
Num_pulse = 30
t = [0,30]  # t(時間)が0〜100まで動き、その時のfを求める。
t_eval = np.linspace(*t, 300)
for i in range(Num_pulse):
    tempul = [i,i+0.8,[0,-1,0]]
    pulse.append(tempul)

print(pulse[0][0])

print(t_eval)
new_pulse = []

for t in t_eval:
    flag = 1
    for Bi in pulse:
        if Bi[0] <= t and Bi[1] >= t:
            new_pulse.append([t,Bi[2]])
            flag = 0
            break
    if flag:
        new_pulse.append([t, [0,0,0]])

print(new_pulse)

