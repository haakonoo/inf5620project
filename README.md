inf5620project
==============

I have based my project on the code in nsbench found on 'http://launchpad.net/nsbench/'
with many modifications.

To run program use ./ns
Specify what kind of parameters you want by comandline arguments.
Example:
./ns womersley2d ipcs dt=0.01 N=100 save_solution=TRUE save_POD=5

Solution will be saved in /results aswell as the 5'th order POD.

For the cylinders N=i implie the mesh cylinder_ik.xml.gz found in /mesh
Complete list of arguments:
"N": 4,
"dt": 0.1,
"save_solution": False,
"save_frequency": 1,
"check_frequency": 10,
"save_solution_at_t=T": False,
"save_xml": False,
"compute_stress": False,
"compute_divergence": False,
"compute_POD": 0,
"compare_POD": 0,
"save_POD": False,
"krylov_solver_absolute_tolerance": 1e-25,
"krylov_solver_relative_tolerance": 1e-12,
"krylov_solver_monitor_convergence": False
