refjast="reference.jast"


##################################################
##GAMESS INPUT
#################################################
gamess_header="""\
 $CONTRL SCFTYP=UHF MULT=1 ICHARG=0 RUNTYP=GRADIENT MAXIT=70 ECP=READ $END
 $CONTRL UNITS=BOHR $END
 $SYSTEM KDIAG=3 $END
 $SCF DIRSCF=.T. RSTRCT=.T. SOSCF=.F.  $END
"""

gamess_guess_line=" $GUESS GUESS=MOREAD NORB=18 $END\n"

gamess_symm="C1"

##################################################

basis=["""\
S 7 1.
   1       25.967331          3.54063093E-03
   2       9.1676598          1.67008793E-02
   3       2.8460805          6.08226508E-02
   4       1.0048708          1.72330686E-01
   5       0.38502821         3.53696261E-01
   6       0.16005760         3.89585261E-01
   7       0.69101229E-01     1.46735635E-01
S   1
  1       0.3258000000          1.000000000
S   1
  1       0.1027000000          1.000000000
P 4 1.
  1      18.1264860      0.002058194485
  2       4.4559865      0.015034403850
  3       1.1582862      0.151153254935
  4       0.3505376      0.890125250947
P 1 1.
  1       .18            1.
"""]

##################################################

ecp_sec="""\
 $ECP
H-TR4 GEN 0   1
6    P
 -9.39168039 2  27.916162
 -0.536989635 2  394.363177
  0.717083592 2  15.663672
 -1.25803771 2  90.1659831
 -0.267890288 2  2232.29185
  1. 1  21.2435951
1    S
   0.        2   1.
H-TR4
 $END
"""

##################################################

class Atom:
    basis=0
    label=""
    pos=[0, 0,0]
    velocity=[0,0,0]
    atomnum=0
    atomweight=0
    def __init__(self, lab, bas, atnum, atweight):
        self.basis=bas
        self.label=lab
        self.atomnum=atnum
        self.atomweight=atweight

##################################################
        
class gos_static_options:
    #general
    eref=-1.1
    vmc_timestep=1.0
    vmc_ndecorr=2
    vmc_nstep=10
    ref_config="reference.config"

    #For optimization
    nconfig=100
    eref=-1.1
    iterations=50
    ref_jast="reference.jast"

    md_timestep=.1
    md_nsteps=1
    md_vmc_nstep=4000
    md_vmc_nblock=20
    md_vmc_nconfig=1
    
static_options = gos_static_options()
    
    

##################################

########################################################
# Start program here
########################################################
#Define the starting position here
atoms=[Atom("H",0,1,1), Atom("H", 0,1,1)]
atoms[0].pos=[0,0,0]
atoms[1].pos=[0,0,1.5]
simulation_time=2

gamess_base="start"
gos_base="startq"
log_file="md.log"

