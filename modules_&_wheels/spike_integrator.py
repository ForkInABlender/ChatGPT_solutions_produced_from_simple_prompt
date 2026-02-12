# Dylan Kenneth Eliot

"""

This file has been moved here for ease of access and use in mapping neural spikes based on mol weight and crippen values plus gate voltage.
"""

class Integrator:
    def __init__(self, y0, deriv, t, dt):
        self.deriv = deriv
        self.y = y0[:]  # Single neuron's state
        self.dt = dt
        self.t = t
        # Initialize with random values between -1.0 and 1.0
        self._w = [random.randrange(-1000000000000000, 1000000000000000)/1000000000000000 for _ in range(len(y0))]
        self._k1 = [random.randrange(-1000000000000000, 1000000000000000)/1000000000000000 for _ in range(len(y0))]
        self._k2 = [random.randrange(-1000000000000000, 1000000000000000)/1000000000000000 for _ in range(len(y0))]
        self._k3 = [random.randrange(-1000000000000000, 1000000000000000)/1000000000000000 for _ in range(len(y0))]
        self._k4 = [random.randrange(-1000000000000000, 1000000000000000)/1000000000000000 for _ in range(len(y0))]

    def step(self):
        self.deriv(self._k1, self.y, self.t)
        for i in range(len(self.y)):
            self._w[i] = self.y[i] + 0.5 * self.dt * self._k1[i]
        self.deriv(self._k2, self._w, self.t + 0.5 * self.dt)
        for i in range(len(self.y)):
            self._w[i] = self.y[i] + 0.5 * self.dt * self._k2[i]
        self.deriv(self._k3, self._w, self.t + 0.5 * self.dt)
        for i in range(len(self.y)):
            self._w[i] = self.y[i] + self.dt * self._k3[i]
        self.deriv(self._k4, self._w, self.t + self.dt)
        dto6 = self.dt / 6.0
        for i in range(len(self.y)):
            self.y[i] += dto6 * (self._k1[i] + 2*self._k2[i] + 2*self._k3[i] + self._k4[i])
        self.t += self.dt
        return self

    def steps(self, n):
        for _ in range(n):
            self.step()
        return self

def IntegratorFactory(y0, deriv, t, dt):
    return Integrator(y0, deriv, t, dt)


"""
Below is example usage. 

def example_deriv(out, y, t):
    for i in range(len(y)):
        out[i] = (molwt + logp - y[i]) / tau # "dv/dt = (molwt + logp - v) / tau" 

molwt= ``['amw']``             # molecular weight term
logp=``['CrippenClogP']``     # partition coefficient
tau     = 1.0              # time constant
v_th    = 0.8              # threshold voltage (Brian2-style)
v_reset = 0.0              # reset voltage
n_neurons = 2              # Number of neurons (to track spikes for each one)
y0 = [random.randrange(-1000000000000000, 1000000000000000)/1000000000000000 for _ in range(n_neurons)] 
t0 = 0.0
dt = 0.01
n_steps = 1000  # simulate for 10 ms if dt=0.01
spikes=[[] for _ in range(n_neurons)]
neurons=[[] for _ in range(n_neurons)]
for _ in range(n_steps):
    for neuron_idx in range(n_neurons):
        integrator = IntegratorFactory([random.randrange(-1000000000000000, 1000000000000000)/1000000000000000], example_deriv, t0, dt)
        integrator.step()
        v = integrator.y[0]
        t = integrator.t        
        if v >= v_th:
            spikes[neuron_idx].append(v)
            integrator.y[0] = v_reset  # Reset voltage
        neurons[neuron_idx].append(integrator)
"""
