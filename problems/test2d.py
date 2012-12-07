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
    "System test"
    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Set parameters
        self.R = 1.21
        self.L = 4
        self.A = 1.3
        self.a = -0.1
        self.n = 2

        # Make channel mesh
        # Divide N by two since the number of triangles will be of order 2*N
        N = int(options["N"]/2)
        self.mesh = Rectangle(0,-self.R,self.L,self.R,int(N*self.L),int(N*2*self.R))

        # Mark bondary
        self.boundary_markers = MeshFunction("uint",self.mesh,1)
        self.boundary_markers.set_all(0)
        NoSlipBoundary(self.R).mark(self.boundary_markers,1)
        InflowBoundary().mark(self.boundary_markers,2)
        OutflowBoundary(self.L).mark(self.boundary_markers,3)

        # Set viscosity (Re = 8)
        self.nu = 1.0 / 8.0

        # Set force
        dpdx = "-sin(8*pi*t)"
        fx = "(a+nu*pi*pi*n*n/(R*R))*A*exp(a*t)*sin(pi*n*x[1]/R)" + dpdx
        self.f =  Expression((fx,"0"), R = self.R, A = self.A, a = self.a, 
                n = self.n, nu = self.nu, t=0)

        # Set exact expression
        self.u_exact = Expression(("A*exp(a*t)*sin(pi*n*x[1]/R)","0"), 
                R = self.R, A = self.A, a = self.a, n = self.n, t=0)
        
        # Set end-time
        self.T = 1

    def initial_conditions(self, V, Q):

        u0 = self.u_exact
        p0 = self.pressure(0,Q)

        return u0, p0

    def boundary_conditions(self, V, Q, t):
        bc = self.boundary_markers

        # Create no-slip boundary condition for velocity
        bv0 = DirichletBC(V, Constant((0.0, 0.0)), bc, 1)
        
        # Create boundary conditions for pressure
        bp0 = DirichletBC(Q, self.pressure(t,Q), bc, 2)
        bp1 = DirichletBC(Q, self.pressure(t,Q), bc, 3)

        bcu   = [bv0]
        bcp   = [bp0,bp1]

        return bcu, bcp

    def pressure(self, t, Q):
        element = FiniteElement("CG", triangle, 1)
        a = -1;
        b = 0;
        return Expression("a*sin(8*pi*t)*x[0] + b",t = t, a = a, b=b, element=element)

    def update(self, t, u, p):
        self.f.t = t

    def __str__(self):
        return "Channel"
