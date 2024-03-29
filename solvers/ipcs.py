# coding: latin1
__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.
# Modified by Håkon Østerbø, 2012

from solverbase import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FunctionSpace(mesh, "DG", 0)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = interpolate(u0, V)
        u1 = Function(V)
        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        problem.f.t = dt
        f  = problem.f
        n  = FacetNormal(mesh)

        # Tentative velocity step
        F1 = inner(u,v)*dx + k*inner(grad(u)*u0,v)*dx \
            + k*2*nu*inner(epsilon(u),grad(v))*dx - k*nu*inner(grad(u).T*n,v)*ds \
            -k*inner(p0,div(v))*dx + k*inner(p0*v,n)*ds - inner(u0,v)*dx - k*inner(f,v)*dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Pressure correction
        a2 = inner(grad(q), grad(p))*dx
        L2 = inner(grad(q), grad(p0))*dx - (1/k)*q*div(u1)*dx

        # Velocity correction
        a3 = inner(v, u)*dx
        L3 = inner(v, u1)*dx - k*inner(v, grad(p1 - p0))*dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3) 

        # Time loop
        self.start_timing()
        for t in t_range:

            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Compute tentative velocity step
            b = assemble(L1)
            [bc.apply(A1, b) for bc in bcu]
            solve(A1, u1.vector(), b, "gmres", "ilu")

            # Pressure correction
            b = assemble(L2)
            [bc.apply(A2, b) for bc in bcp]
            solve(A2, p1.vector(), b, 'gmres', 'hypre_amg')

            # Velocity correction
            b = assemble(L3)
            [bc.apply(A3, b) for bc in bcu]
            solve(A3, u1.vector(), b, "gmres", "ilu")

            # Update
            self.update(problem, t, u1, p1)
            u0.assign(u1)
            p0.assign(p1)
        return u1, p1

    def __str__(self):
        return "IPCS"
