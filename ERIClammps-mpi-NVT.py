#!/usr/bin/python
#
# Python script that runs a LAMMPS simulation.  Yay, the best of both worlds!
#
# Generally, you want to run this script in parallel.  In that case, first run
# in bash:
#
#     module add openmpi-x86_64
#
# and then execute:
#
#     mpirun  -np $NP  --map-by core  lammps-run-NVT.py  arg1  arg2  ...
#
# where $NP is the number of cores used, "lammps-run-NVT.py" the name of this
# script, and "arg1", "arg2", etc. the arguments of this script.
#
#
#                                      by Erik Lascaris, version 27 June 2016


# -----------------------------------------------------------------------------

from math import *
from sys import argv, exit
from time import strftime  # for timestamp


# ----  argument parser  ------------------------------------------------------

import argparse

# Read the arguments.
parser = argparse.ArgumentParser(description='Run a standard LAMMPS simulation (NVT).')

"""
parser.add_argument('-f1', dest='file1', action='store', nargs=1, required=True, #type=file,
                    help='1st configuration file (*.gro, *.g96, etc.)')

parser.add_argument('-dr1', dest='dr1', action='store', nargs=3, type=float, default=[0,0,0],
                    help='move all atoms in 1st config file by  dx dy dz  (values in nm)')


parser.add_argument('-f2', dest='file2', action='store', nargs=1, required=True, #type=file,
                    help='2nd configuration file (*.gro, *.g96, etc.)')

parser.add_argument('-dr2', dest='dr2', action='store', nargs=3, type=float, default=[0,0,0],
                    help='move all atoms in 2nd config file by  dx dy dz  (values in nm)')


parser.add_argument('-o', dest='out', action='store', nargs=1, required=True, #type=argparse.FileType('w'),
                    help='name of output configuration file (*.gro, *.g96, etc.)')


parser.add_argument('-box', dest='box', action='store', nargs=3, type=float, default=[0,0,0],
                    help='set box of output file to  Lx Ly Lz  (values in nm)')

parser.add_argument('-p', dest='topol', action='store', nargs=1, default="",
                    help='name of topol file that will be updated (typically "topol.top")')


args = parser.parse_args()
"""
#print args
#print


#    fout = args.out[0]




# ----  start MPI code  -------------------------------------------------------

from mpi4py import MPI

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

master = (rank == 0)

if (master):
    #sys.stdout.write(
    print "Running LAMMPS using mpi4py !!  This is process %d of %d on %s.\n" % (rank, size, name) #)





# TODO --- catch error somehow if we don't have MPI available (b/c of module for example)
#import pypar

# See documentation at:
#   http://www.shocksolution.com/microfluidics-and-biotechnology/parallel-programming-with-python/pypar-documentation/

#
#if (pypar.initialized()):
#    rank = pypar.rank()
#if (rank == 0): print 'MPI is working!!! :-)  Using {0} processors'.format(pypar.size())
#mpi_initialized = True

#else:
#    rank = 0
#    print 'possible issue with MPI...'
#    print "If running parallel, do this first:  module add openmpi-x86_64"
#    pypar_initialized = False


#master = (rank == 0)


from lammps import lammps

lmp = lammps()



if (master):
    #print 'run as:    mpirun -np 48 '+argv[0]+'  crystal.cfg  2000'
    #print 'run as:    mpirun -np 48 '+argv[0]+'  init-config.cfg  2000  [opt:runtime-ps]'
    print "Run as:   mpiexec  -n 48  python {0}  init.config.cfg  2000  [opt:runtime-ps]".format(argv[0])


# Check number of args
if (len(argv)-1 < 2):
    if (master):
        print
        print 'ERROR! Need 2 arguments: input file (cfg), and temperature (in K).'
        print 'Optional: 3rd argument with run time in ps.'
        print
    #if (pypar_initialized):  pypar.finalize()
    exit(1)



fname_in = argv[1]
T = float(argv[2])



if (len(argv) >= 3+1):
    run_time = int(1000 * float(argv[3]) + 0.5)  # +0.5 for rounding
else:
    run_time = 10000    # good test run; 10 ps
    #run_time = 100000   # minimum production run; 100 ps
    #run_time = 1000000  # serieus run; 1 ns



if (master):
    #print ''
    print 'input file  =  ' + fname_in
    print 'temperature =  ' + str(T) + ' K'
    print "Run time    =  {0} ps".format(run_time/1000)

#run_time = 600000  # XXX 10000 = 10 min [so 10 hours = 600 min => 600,000]

# XXX -- set at top of script

#run_time = 10000 #00  # NEEDS TO BE 100 t0 = ABOUT 100 ns = 100,000,000   /  150 sec for 10,000  so 3 hours = 720,000 --> 1,000,000 (4 hours) on 48 cores

nsteps_out = run_time/2000  # 500  Note: 10000 = is ok

#nsteps_traj = run_time/200  #20000    #20   # NOTE: 20000 => 51GB, so 200 => 510 MB.

# for 3456 atoms (12^3 box) run_time/200 gives 26M traj file.
nsteps_traj = 100
#nsteps_traj = run_time/10000
#nsteps_traj = run_time/200    



# Ensure proper values, of at least 1.
if (nsteps_out == 0):   nsteps_out = 1
if (nsteps_traj == 0):  nsteps_traj = 1


# 

#lmp.file("input.inp")  # TODO -- do all LAMMPS commands in this file, and load here instead a configuration file.
#lmp.file('Relax.inp')



# Initialize the LAMMPS simulation.
lmp.command('units       metal')
lmp.command('dimension   3')
lmp.command('boundary    p  p  p')
lmp.command('newton      on')
lmp.command('atom_style  atomic')


# IF WE NEED TO READ IN A CONFIG FILE, DO THIS:

#if (False        and  read_config_file):  # TODO

lmp.command('read_data {0}'.format(fname_in)) #crystal.cfg')




# Tell LAMMPS what force field to use.
lmp.command('pair_style eam/cd')
lmp.command('pair_coeff * * FeCr.cdeam Fe Cr')
lmp.command('mass 1 55.8470')   # type 1 = Fe  (atom nr 26)
lmp.command('mass 2 51.9961')   # type 2 = Cr  (atom nr 24)

#
lmp.command('neighbor 1.0 bin')
lmp.command('neigh_modify every 1 delay 10 check yes')

# T = 100 K!!
lmp.command('velocity all create {0} 87289 rot yes dist gaussian'.format(T))


# temp values = Tstart Tstop Tdamp (typical: Tdamp = 100 timesteps = 0.1 ps)
# iso values = Pstart Pstop Pdamp  (typical: Pdamp = 1000 timesteps = 1 ps)
#lmp.command('fix 1 all   npt   temp  {0} {0}  0.1    iso  0.0  0.0  0.1'.format(T))
lmp.command('fix 1 all   nvt   temp  {0} {0}  0.1'.format(T))

#lmp.command('compute Ekin all ke/atom')   # to keep track of fastest atom, so we can adjust dt.')
#lmp.command('compute maxEkin all reduce max c_Ekin')


lmp.command('group GRP_Fe type 1')
lmp.command('compute msdFe GRP_Fe msd')

lmp.command('group GRP_Cr type 2')
lmp.command('compute msdCr GRP_Cr msd')

lmp.command('group GRP_All subtract GRP_Cr GRP_Fe')

#lmp.command('compute cnaAll all cna/atom 1.0') #3.45')
#lmp.command('compute centro all centro/atom bcc')


lmp.command('thermo {0}'.format(nsteps_out))  # 100 = how often you output')  -- XXX  changed from 100 to 200.
lmp.command('thermo_style custom  time temp press vol pe ke enthalpy c_msdFe[4] c_msdCr[4]')
lmp.command('thermo_modify norm no flush yes')


# XXX --- problem with VMD reading XTC [you need GRO file too..]
#lmp.command('dump DUMP1 all xtc 2000 traj.xtc')  # XXX at 200 (and run=10000) we get 27 MB.  So we'll do 2000 which reduces 10x size [so run=600,000 would give 150 MB]
lmp.command('dump DUMP1 all custom {0}  traj.lammpstrj  element id x y z'.format(nsteps_traj))
#lmp.command('dump DUMP1 all custom {0}  traj.lammpstrj  element id x y z c_cnaAll c_centro'.format(nsteps_traj))  # when 200 then 152 MB for run=6000.. so we'll do 100x
lmp.command('dump_modify DUMP1 element Fe Cr')
lmp.command('dump_modify DUMP1 sort id')


# Set the timestep.  This is also used to determine the status while running.
lmp.command('timestep 0.001')  # 0.001 ps = default for units=metal

# Before the run, write a timestamp to the log file, so that get_status.sh can
# determine when the run started (and estimate when it will finish).
lmp.command('print "TIME= {0}"'.format(strftime("%Y-%m-%d %H:%M:%S")))

# Run the LAMMPS simulation.
lmp.command('run {0} pre yes post yes'.format(run_time))

# Write final configuration file.
lmp.command('write_data  out.cfg')


# Probably not needed because script ends, but just to be proper:
#if (pypar_initialized):
#    pypar.finalize()
#    quit()




