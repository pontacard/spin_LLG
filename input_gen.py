import numpy as np
from scipy import signal

class InputGenerator:

    def __init__(self, start_time, end_time, num_time_steps):
        self.start_time = start_time
        self.end_time = end_time
        self.num_time_steps = num_time_steps

    def generate_sin(self, amplitude=1.0):
        t = np.linspace(self.start_time, self.end_time, self.num_time_steps, endpoint=False)
        y = -(signal.square(2 * np.pi *0.1 * t)+1)/2
        return y
