class Rossler:
    '''
        Initializes a Rössler-System described by the following system of differential equations:
        dx/dt = -y - z
        dy/dt = x + ay
        dz/dt = b + z(x - c)

        Attributes
        __________

        x0, y0, z0 : Floats
            The initial values.
        a, b, c : Floats
            The parameters of the dynamical system. Different parameters produce different dynamics.
            The default parameters produce the famous strange attractor.

        Methods
        _______

        dx, dy, dz
            Implementation of the differential equations. Describe the change of the system, given a state.
        integrate
            Evolves the system with the above differential equations, using a 4th order Runge-Kutta scheme.
    '''
    def __init__(self, x0=1, y0=1, z0=0, a=0.2, b=0.2, c=5.7):
        self.a = a
        self.b = b
        self.c = c
        self.x = x0
        self.y = y0
        self.z = z0
        self.trajectory = [tuple([x0, y0, z0])]

    # The differential equations
    def dx(self, x, y, z):
        return -y - z

    def dy(self, x, y, z):
        return x + self.a * y

    def dz(self, x, y, z):
        return self.b + z * (x - self.c)

    '''
        Implementation of a 4th order Runge-Kutta integration scheme. Takes the current state of the system and evolves
        it according to the differential equations (defined above as dx, dy, dz).
        
        Parameters
        __________
        
        steps : Int
            The number of times the integration procedure is applied to the system.        
        h : Float
            Integration step size. The bigger h, the bigger the error.
    '''
    def integrate(self, steps, h=0.01):
        for step in range(steps):
            k1x = h * self.dx(self.x, self.y, self.z)
            k1y = h * self.dy(self.x, self.y, self.z)
            k1z = h * self.dz(self.x, self.y, self.z)

            k2x = h * self.dx(self.x + 0.5 * k1x, self.y + 0.5 * k1y, self.z + 0.5 * k1z)
            k2y = h * self.dy(self.x + 0.5 * k1x, self.y + 0.5 * k1y, self.z + 0.5 * k1z)
            k2z = h * self.dz(self.x + 0.5 * k1x, self.y + 0.5 * k1y, self.z + 0.5 * k1z)

            k3x = h * self.dx(self.x + 0.5 * k2x, self.y + 0.5 * k2y, self.z + 0.5 * k2z)
            k3y = h * self.dy(self.x + 0.5 * k2x, self.y + 0.5 * k2y, self.z + 0.5 * k2z)
            k3z = h * self.dz(self.x + 0.5 * k2x, self.y + 0.5 * k2y, self.z + 0.5 * k2z)

            k4x = h * self.dx(self.x + k3x, self.y + k3y, self.z + k3z)
            k4y = h * self.dy(self.x + k3x, self.y + k3y, self.z + k3z)
            k4z = h * self.dz(self.x + k3x, self.y + k3y, self.z + k3z)

            xx = self.x + (1 / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
            yy = self.y + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)
            zz = self.z + (1 / 6) * (k1z + 2 * k2z + 2 * k3z + k4z)

            self.x = xx
            self.y = yy
            self.z = zz

            self.trajectory.append(tuple([xx, yy, zz]))

    def clear(self):
        self.trajectory = []

    def get_x(self):
        return list(zip(*self.trajectory))[0]

    def get_y(self):
        return list(zip(*self.trajectory))[1]

    def get_z(self):
        return list(zip(*self.trajectory))[2]

