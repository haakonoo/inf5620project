__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2010-04-15"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# List of solvers
solvers = ["ipcs"]

# Wrapper for solver classes
def Solver(name, options):
    "Return solver instance for given solver name"
    exec("from %s import Solver as NamedSolver" % name)
    return NamedSolver(options)
