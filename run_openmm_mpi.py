from mpi4py import MPI
import os
import subprocess
import argparse

''' 
Start mpi process
'''

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

'''Arg parsing'''

parser = argparse.ArgumentParser()
parser.add_argument('-T',
                    '--type_sim',
                    type=str,
                    help='Type of simulation: eq, prod0, prod1...')

parser.add_argument('-P',
                    '--PWD',
                    type=str,
                    help='working directory')

parser.add_argument('-S',
                    '--structdir',
                    type=str,
                    help='directory with structural files')

parser.add_argument('-o',
                    '--outdir_gen',
                    type=str,
                    help='general directory for output replica')

parser.add_argument('-n',
                    '--numsteps',
                    type=int,
                    required=False,
                    default=100000000,
                    help='Number of production steps')

args = parser.parse_args()

if rank<4:
    try:
        os.mkdir(f'{args.PWD}/{args.outdir_gen}')
    except:
        pass

comm.barrier()


try:
    os.mkdir(f'{args.PWD}/{args.outdir_gen}/replica{rank}')
except:
    print(f"replica{rank} already exists probably")



''' Depending on type of simulation run different jobs 
    eq: equlibration
    prod0: production 0
    prod1: production 1
    ...
    write a new if statement 
        to run another restart siulation

    Future: set up code to take in multiple structure files 
'''

if args.type_sim=='eq':
    from src.simulate.equilibrate import *

    inpcrd_fil = f'{args.PWD}/{args.structdir}/inpcrd1'
    prmtop_fil = f'{args.PWD}/{args.structdir}/prmtop1'
    try:
        os.mkdir(f'{args.PWD}/{args.outdir_gen}/replica{rank}/eq')
    except:
        pass
                                           
    eq_st = f'{args.PWD}/{args.outdir_gen}/replica{rank}/eq/eq.state'
    eq_chkpt = f'{args.PWD}/{args.outdir_gen}/replica{rank}/eq/eq.chk'
    device = f'{rank%4+2}'
    eq_simulation = running_eq(inpcrd_fil,
                                prmtop_fil,
                                eq_st,
                                eq_chkpt,
                                device)

elif args.type_sim=='prod0':
    from src.simulate.production_0 import *
    try:                                              
        os.mkdir(f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0')
    except:
        pass

    inpcrd_fil = f'{args.PWD}/{args.structdir}/inpcrd1'
    prmtop_fil = f'{args.PWD}/{args.structdir}/prmtop1'
    device = f'{rank%4+2}'
    eq_chkpt = f'{args.PWD}/{args.outdir_gen}/replica{rank}/eq/eq.chk'
    prod_steps = args.numsteps #100000000
    prod_dcd = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0/prod0.dcd'
    prod_rst = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0/prod0.rst.chk'
    prod_chkpt = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0/prod0.chk'
    prod_st = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0/prod0.state'
    prod_log = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0/prod0.csv'
    
    prod0(inpcrd_fil,
      prmtop_fil,
      device,
      eq_chkpt,
      prod_steps,
      prod_dcd,
      prod_rst,
      prod_chkpt,
      prod_st,
      prod_log)

elif args.type_sim=='prod1':
    from src.simulate.production_1 import *
    try:                                              
        os.mkdir(f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod1')
    except:
        pass

    inpcrd_fil = f'{args.PWD}/{args.structdir}/inpcrd1'
    prmtop_fil = f'{args.PWD}/{args.structdir}/prmtop1'
    device = f'{rank%4+2}'
    eq_chkpt = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod0/prod0.rst.chk'
    prod_steps = args.numsteps #100000000
    prod_dcd = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod1/prod1.dcd'
    prod_rst = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod1/prod1.rst.chk'
    prod_chkpt = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod1/prod1.chk'
    prod_st = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod1/prod1.state'
    prod_log = f'{args.PWD}/{args.outdir_gen}/replica{rank}/prod1/prod1.csv'
    
    prod1(inpcrd_fil,
      prmtop_fil,
      device,
      eq_chkpt,
      prod_steps,
      prod_dcd,
      prod_rst,
      prod_chkpt,
      prod_st,
      prod_log)





    




if False:
    # Print the current working directory
    #print("Current Directory:", os.getcwd())
    #PWD = '/eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_Dimer'
    #replicas = [rank*4+it for it in range(1,5)]
    for r, d in zip(replicas,[0,1,2,3]):
        # Change the current working directory
        try:
            os.mkdir(f"{PWD}/replica{r}/")
    
        except:
            pass
    
        os.chdir(f"{PWD}/replica{r}")
        print("Directory after change:", os.getcwd())
    
        # Define the CPU cores you want to pin the process to (e.g., cores 0 and 1)
        #cpu_cores = [d*8+it for it in range(8)]
        
        # Set the CPU affinity
        #os.sched_setaffinity(d, cpu_cores)
        
        process = subprocess.Popen(["python", 
                                    "equilibrate.py",
                                    "-R",
                                    f"{PWD}/replica{r}",
    
                                    ">", "equilibrate.log", "2>", "equilibrate.err", "&"])
        processes.append(process)
    #    else:
    #        process = subprocess.Popen(["python", "production.1.py", ">", "production.1.log", "2>", "production.1.err", "&"])
    #        processes.append(process)
    #
    #for process in processes:
    #    process.wait()
    
    print("All scripts have finished.")
