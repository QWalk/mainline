import os
import string

gamess2qmc="gamess2qmc"
gamess="nohup nice /home/apps/gamess/rungms"
gos="/home/lkwagner/qmcsources/qmc_inherit/gos-Linux-3.2"


#We want to read the sysdef module from the current
#directory, so add this to Python's search path
import sys
sys.path.insert(0,".")

from sysdef import *

def print_gamess(atoms, file, basis, ecp_section, last_pun):
    f=open(file, "w")
    f.write(gamess_header)
    if(last_pun != ''):
        f.write(gamess_guess_line)
    f.write(" $DATA \n \n")
    f.write(gamess_symm+"\n")
    for at in atoms:
        string=" "+at.label+" "+ str(at.atomnum) + "  " \
                + str(at.pos[0]) + "  " + str(at.pos[1]) +  "  "\
                + str(at.pos[2]) + "\n"
        f.write(string);
        f.write(basis[at.basis]+ "\n")
    f.write(" $END\n")
    f.write(ecp_section)
    if(last_pun != ''):
        inf=open(last_pun, "r")
        pr=0; 
        while(1):
            line=inf.readline()
            if(line==''):
                break
            if(line.rfind("$VEC") > 0) or line.rfind("$SCF") > 0:
                pr=1
            if(pr==1):
                f.write(line)
    f.close()

##################################################

class gos_dynamic_options:
    root=""
    wffile=""
    md_wf=""
    md_read=""
    md_write=""


###################################################

def print_gos_opt(sopts, dopts):
    "Print optimization file"
    dopts.md_wf=dopts.root+".wf"
    return """\
METHOD {
  VMC
  NBLOCK 1
  NSTEP 5
  TIMESTEP """ + str(sopts.vmc_timestep) + """
  NDECORR """ + str(sopts.vmc_ndecorr) +  """
  STORECONFIG md_opt.config
  READCONFIG md_opt.config
  NCONFIG """ + str(sopts.nconfig) + """
  EREF """ + str(sopts.eref) + """
}

METHOD {
  OPTIMIZE
  ITERATIONS """ + str(sopts.iterations)+ """
  READCONFIG  md_opt.config
  NCONFIG """ + str(sopts.nconfig) + """
  EREF """ + str(sopts.eref) + """
  MINFUNCTION VARIANCE
  WFOUTPUT """ + dopts.md_wf + """
  PSEUDOTEMP /scratch/lkwagner/qmcpseudo
}

INCLUDE """ + dopts.root + """.sys
TRIALFUNC { INCLUDE """ + dopts.wffile + " } \n"   


##################################################

def print_gos_md(sopts, dopts):
    read=''
    write=''
    if(dopts.md_read != ''):
        read="READCHECK " + dopts.md_read + "\n"
    if(dopts.md_write !=''):
        write="STORECHECK " + dopts.md_write +"\n"
    atomweights=''
    for at in atoms:
        atomweights=atomweights + str(at.atomweight) + " "
            
    return """\
METHOD {
  MD
  NSTEP """ +str(sopts.md_nsteps)+ """
  TIMESTEP """ + str(sopts.md_timestep) + "\n  " + \
   read + write +  """
  ATOMIC_WEIGHTS { """ + atomweights + """ } 
  EMBED {
   VMC
   NBLOCK   """ + str(sopts.md_vmc_nblock) + """
   NSTEP    """ + str(sopts.md_vmc_nstep) + """
   NDECORR  """ + str(sopts.vmc_ndecorr) + """
   TIMESTEP """ + str(sopts.vmc_timestep) + """
   NCONFIG  """ + str(sopts.md_vmc_nconfig) + """
   STORECONFIG  md.config 
   """ + sopts.md_vmc_extra + """
   #uncomment the following to read a configuration
   READCONFIG  md.config
   EREF """ + str(sopts.eref) + """
  }
}

INCLUDE """ + dopts.root + """.sys
TRIALFUNC { INCLUDE """ + dopts.md_wf+ " } \n"    


##################################################
##################################################



#log=open(log_file, "w")

nsteps=int(simulation_time/(static_options.md_timestep*static_options.md_nsteps))

for i in range(0,nsteps):

    #Run GAMESS
    gamessroot=gamess_base + string.zfill(str(i),4)
    gosroot=gos_base+string.zfill(str(i),4)
    last_pun=''
    if(i > 0):
        last_pun=gamess_base+string.zfill(str(i-1), 4) + ".pun"
    elif(gamess_punch_start!= ''):
        last_pun=gamess_punch_start
        
    print_gamess(atoms, gamessroot + ".inp", basis, ecp_sec, last_pun)
    systext=gamess+" " + gamessroot + ' >& ' + gamessroot + ".out"
    print "executing: " + systext
    os.system(systext)

    #Convert from GAMESS to QMC
    
    systext=gamess2qmc +" -o " + gosroot + " " + gamessroot
    print "executing: " + systext
    os.system(systext)
    os.remove(gosroot+".hf");
    os.remove(gosroot+".jast");
    os.remove(gosroot+".jast2");
    
    refwf=gosroot+".refwf"
    
    f=open(refwf, "w");
    f.write("SLATER-JASTROW\n  WF1 { INCLUDE " +\
            gosroot + ".slater }\n  WF2 { INCLUDE " + refjast + " }\n")
    f.close()
    
    
    #Set dynamic options

    firstoptfile=gosroot+".opt"
    dyn_opt=gos_dynamic_options()
    dyn_opt.root=gosroot
    dyn_opt.wffile=refwf
    dyn_opt.md_write=gosroot+".mdcheck"
    if(i > 0):
        last_name=gos_base+string.zfill(str(i-1),4)
        dyn_opt.md_read=last_name + ".mdcheck"
    

    #Optimize the first wave function

    f=open(firstoptfile, "w")
    f.write(print_gos_opt(static_options, dyn_opt))
    f.close()
    
    systext=gos + " " + firstoptfile
    print "executing: " + systext
    os.system(systext)

    #Run MD
    
    mdfile=gosroot+".md"
    f=open(mdfile, "w")
    f.write(print_gos_md(static_options, dyn_opt))
    f.close()
    
    systext=gos + " " + mdfile
    print "executing: " + systext
    os.system(systext)

    

    mdout=open(dyn_opt.md_write, "r")
    mdout.readline() #natoms
    mdout.readline() #current_pos
    n=0
    for at in atoms:
        line=mdout.readline()
        lsp=line.split()
        at.pos[0]=float(lsp[0])
        at.pos[1]=float(lsp[1])
        at.pos[2]=float(lsp[2])

