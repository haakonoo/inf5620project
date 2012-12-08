# coding: latin1
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-04-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Håkon Østerbø, 2012

from dolfin import *

import numpy as np

from time import time
from os import getpid
from commands import getoutput

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4

class SolverBase:
    "Base class for all solvers."

    def __init__(self, options):

        # Store options
        self.options = options

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        # Reset files for storing solution
        self._ufile = None
        self._pfile = None

        # Reset storage for functional values and errors
        self._t = []
        self._M = []
        self._m = []
        self._e = []

        # Store solution for POD computation
        self.us = []

        # Store POD basis and coefficiants
        self.A = []
        self.ub = []



    def getMyMemoryUsage(self):
        mypid = getpid()
        mymemory = getoutput("ps -o rss %s" % mypid).split()[1]
        return mymemory

    def start_timing(self):
        """Start timing, will be paused automatically during update
        and stopped when the end-time is reached."""
        self._time = time()

    def solve(self, problem, dt, plot_solution=True):
        "Solve problem"
        raise NotImplementedError

    def prefix(self, problem):
        "Return file prefix for output files"
        p = problem.__module__.split(".")[-1].lower()
        s = self.__module__.split(".")[-1].lower()
        return problem.output_location + p + "_" + s

    def update(self, problem, t, u, p):
        "Update problem at time t"

        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Compute divergence
        if self.options["compute_divergence"]:
            check_divergence(u, p.function_space())

        # Update problem
        problem.update_problem(t, u, p)

        # Error check
        if problem.u_exact != None:
            problem.u_exact.t = t
            V = u.function_space()
            error = Function(V)
            error.vector()[:] = u.vector() - interpolate(problem.u_exact, V).vector()
            e = norm(error)
        else :
            e = None

        # Store values
        self._t.append(t)
        self._e.append(e)

        # Save solution
        if self.options["save_solution"]:

            # Save velocity and pressure
            frequency = self.options["save_frequency"]
            N = self.options['N']
            if (self._timestep - 1) % frequency == 0:
                # Create files for saving
                if self._ufile is None:
                    self._ufile = File("results/" + self.prefix(problem) +"N_"+ str(N) + "_u.pvd")
                if self._pfile is None:
                    self._pfile = File("results/" + self.prefix(problem) +"N_"+ str(N) + "_p.pvd")
                self._ufile << u
                self._pfile << p

        # Save solution at t = T
        if self.options["save_solution_at_t=T"]:
            if t >= problem.T:
                N = self.options["N"]
                # Create files for saving
                if self._ufile is None:
                    self._ufile = File("results/" + self.prefix(problem) +"_N_"+ str(N) + "_at_end" + "_u.pvd")
                if self._pfile is None:
                    self._pfile = File("results/" + self.prefix(problem) +"_N_"+ str(N) + "_at_end" + "_p.pvd")
                self._ufile << u
                self._pfile << p

        # Save vectors in xml format
        if self.options["save_xml"]:
            file = File(self.prefix(problem) +"_N_"+ str(self.options["N"]) + "t=%1.2e"% t + "_u.xml" )
            file << u.vector()

            file = File(self.prefix(problem) +"_N_"+ str(self.options["N"]) + "t=%1.2e"% t + "_p.xml" )
            file << p.vector()

        if self.options["store_u"]:
            self.us.append(u.copy(True))

        # Print progress
        print ""
        s = "Time step %d finished in %g seconds, %g%% done (t = %g, T = %g)." \
            % (self._timestep, timestep_cputime, 100.0*(t / problem.T), t, problem.T)
        print s + "\n" + len(s)*"-"

        # Increase time step and record current time
        self._timestep += 1
        self._time = time()

    def meanerror(self):
        "Return mean value of error from computation if given exact solution. Else return None"

        # Return value
        if self._e[0] is None:
            return None
        else:
            return 1./len(self._e)*sum(self._e)

    def cputime(self):
        #"Return accumulated CPU time."
        return self._cputime

    def compute_POD(self,k):
        us = self.us 
        V = us[0].function_space()
        M = len(us)

        # Compute mean velocity
        u_0 = Function(V)
        for u in us:
            u_0.vector()[:] += u.vector()
        u_0.vector()[:] = 1./M*u_0.vector()

        # Compute varying velocity
        um = []
        for u in us:
            um_i = u.copy(True)
            um_i.vector()[:] -= u_0.vector()
            um.append(um_i)
        print "Mean and varying velocity computed" 

        # Compute correlation matrix
        C = np.zeros([M,M])
        for n in range(M):
            for m in range(n+1):
                C[n,m] = (1./M)*assemble(inner(um[n],um[m])*dx)
                C[m,n] = C[n,m]
        print "Correlation matrix computed"

        # Compute eigenvalues
        eigval, eigvec = np.linalg.eig(C)
        
        # Remove complex values due to roundoff errors
        eigval = eigval.real
        eigvec = eigvec.real

        # Scale eigenvectors
        A = self.A
        for i in range(M):
            if eigval[i] > 10e-16:
                A.append(np.sqrt(eigval[i]*M)*eigvec[:,i])
            else : 
                break;

        # Assure that order of POD is ok.
        if k > len(A):
            k = len(A)
        print "Eigenvalues/vectors computed and scaled"
        
        # Basis velocity
        ub = self.ub
        for i in range(k):
            phi_i = Function(V)
            for m in range(M):
                phi_i.vector()[:] += A[i][m]*um[m].vector()
            phi_i.vector()[:] = 1./(M*eigval[i])*phi_i.vector()
            ub.append(phi_i)
        print "Basis computed"

        return u_0, k

    def write_POD_velocity_to_file(self,problem, k):
        u_0, k = self.compute_POD(k)
        M = len(self.us)
        N = self.options["N"]
        podfile = File("results/" + self.prefix(problem) +"_POD_order_"+ str(k) + "_N_"+ str(N) + "_u.pvd")
        for m in range(M):
            u_m = u_0.copy(True)
            for i in range(k):
                u_m.vector()[:] += self.A[i][m]*self.ub[i].vector()

            podfile << u_m
        print "POD written to file"

    def POD_error_analysis(self):
        us = self.us
        M = len(us)

        u0, k= self.compute_POD(M)
        E = []
        
        um = [u0.copy(True) for i in range(M)]

        for i in range(k):
            for m in range(M):
                um[m].vector()[:] += self.A[i][m]*self.ub[i].vector()
            E.append(compute_mean_norm(us,um))
        
        return E

def compute_mean_norm(u,v):
    V = u[0].function_space()
    M = len(u)
    mean_norm = 0
    for i in range(M):
        e = Function(V)
        e.vector()[:] = u[i].vector() - v[i].vector()
        mean_norm += norm(e)
    
    return mean_norm/float(M)

def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, nu):
    "Return stress tensor."
    return 2*nu*epsilon(u) - p*Identity(u.cell().d)

