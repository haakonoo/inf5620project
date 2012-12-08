# coding: latin1
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kent-Andre Mardal, 2008.
# Modified by Håkon Østerbø, 2012.

from dolfin import *
from numpy import linspace
from math import *

class ProblemBase:
    "Base class for all problems."

    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh = None
        self.f    = None
        self.bcu  = []
        self.bcp  = []
        self.nu   = None
        self.t    = 0
        self.T    = None
        self.dt   = None
        self.u0   = None
        self.p0   = None
        self.u    = None
        self.p    = None
        self.U    = 1.0
        self.output_location = ''
        self.u_exact = None

        # FIXME: Temporary while testing
        self.bcpsi = []

    def update_problem(self, t, u, p):
        "Update problem at time t"

        # Update state
        self.t = t
        self.u = u
        self.p = p

        # Call problem-specific update
        self.update(t, u, p)

    def update(self, t, u, p):
        "Problem-speficic update at time t"
        pass

    def tolerance(self, problem):
        "Return tolerance (used as local convergence criterion)."
        if str(problem) == 'Channel':
            return 1e-11
        elif str(problem) == 'Cylinder':
            return 1e-7
        else:
            return 1e-6

    def timestep(self, problem):
        "Return time step and number of time steps for problem."

        T = float(problem.T)
        dt = self.options['dt']
        n = int(T/float(dt) +1)

        # Compute range
        t_range = linspace(0,T,n+1)[1:] 

        # Report time step
        print " "
        print 'Number of timesteps:' , len(t_range)
        print 'Size of timestep:' , dt
        print " "

        return dt, t_range[0], t_range
