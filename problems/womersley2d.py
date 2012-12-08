# coding: latin1
__author__ = "Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-04-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Håkon Østerbø, 2012

from problembase import *
from numpy import array

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return (x[0] < 100*DOLFIN_EPS) and on_boundary

# Outflow boundary
class OutflowBoundary(SubDomain):
        def __init__(self,l):
                SubDomain.__init__(self)
                self.l = l
        def inside(self, x, on_boundary):
                return (x[0] > self.l - 100*DOLFIN_EPS) and on_boundary

# No-slip boundary
class NoSlipBoundary(SubDomain):
        def __init__(self,r):
                SubDomain.__init__(self)
                self.r = r;
        def inside(self, x, on_boundary):
                return (x[1] < -self.r + 100*DOLFIN_EPS or x[1] > self.r - 100*DOLFIN_EPS) and on_boundary

# Problem definition
class Problem(ProblemBase):
    "2D channel test problem with known analytical solution."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.l = 4;
        self.r = 1;
        
        # Create mesh
        N = options["N"]
        N = 4
        self.mesh = Rectangle(0,-self.r, self.l, self.r, int(N*self.l),int(N*2*self.r), "crossed")

        # Mark boundary
        self.boundary_markers = MeshFunction("uint",self.mesh,1)
        self.boundary_markers.set_all(0)
        NoSlipBoundary(self.r).mark(self.boundary_markers,1)
        InflowBoundary().mark(self.boundary_markers,2)
        OutflowBoundary(self.l).mark(self.boundary_markers,3)

        # Create right-hand side function with pressure gradient as body force
        self.f = Constant((0, 0))

        # Set viscosity (Re = 8)
        self.nu = 1.0 / 8.0

        self.U = 2.5

        # Set end-time
        self.T = 1

    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = self.pressure(0,Q)

        return u0, p0

    def boundary_conditions(self, V, Q, t):
        bm = self.boundary_markers

        # Create no-slip boundary condition for velocity
        bv0 = DirichletBC(V, Constant((0.0, 0.0)), bm, 1)
        
        # Create boundary conditions for pressure
        bp0 = DirichletBC(Q, self.pressure(t,Q), bm, 2)
        bp1 = DirichletBC(Q, self.pressure(t,Q),  bm, 3)

        bcu   = [bv0]
        bcp   = [bp0,bp1]

        return bcu, bcp

    def pressure(self, t, Q):
        element = FiniteElement("CG", triangle, 1)
        a = -1;
        b = 0;
        return Expression("a*sin(8*pi*t)*x[0] + b",t = t, a = a, b=b, element=element)

    def functional(self, t, u, p):
        if t < self.T:
            return 0
        else:
            x = array((1.0, 0.5))
            values = array((0.0 , 0.0))
            #u.eval(values, x)
            return values[0]

    def reference(self, t):
        if t < self.T:
            return 0
        else:
            num_terms = 10000
            u = 1.0
            c = 1.0
            for n in range(1, 2*num_terms, 2):
                a = 32.0 / (DOLFIN_PI**3*n**3)
                b = (1/8.0)*DOLFIN_PI**2*n**2
                c = -c
                u += a*exp(-b*t)*c
            return u

    def __str__(self):
        return "Channel"
