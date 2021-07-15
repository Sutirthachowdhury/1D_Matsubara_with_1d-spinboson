import random
import shutil
import os


def sbatch(filename):
    """
    Submit a job and get the job id returned
    """
    submit = os.popen("sbatch  %s"%(filename)).read()
    subId = submit.split()[3].replace("\n","")
    return subId

def makeinp(Fold= "./"):
    fob = open(Fold+"/params.in","w+")

    nstep = 3000
    fob.write(str(nstep) +  "\t" + "nstep" + "\n")

    ntraj = 5000
    fob.write(str(ntraj) +  "\t" + "ntraj" + "\n")

    nmc = 1000
    fob.write(str(nmc) +  "\t" + "nmc" + "\n")

    beta = 16.0
    fob.write(str(beta) +  "\t" + "beta" + "\n")

    mnuc = 1.0
    fob.write(str(mnuc) +  "\t" + "mnuc" + "\n")

    k1 = 1.0
    fob.write(str(k1) +  "\t" + "k1" + "\n")

    kappa = 0.1/2
    fob.write(str(kappa) +  "\t" + "kappa" + "\n")

    alpha = -0.5
    fob.write(str(alpha) +  "\t" + "alpha" + "\n")

    delta = 0.10
    fob.write(str(delta) +  "\t" + "delta" + "\n")

    dt = 0.01
    fob.write(str(dt) +  "\t" + "dt" + "\n")

    ourseed = random.randint(0,1E6)
    fob.write(str(ourseed) +  "\t" + "ourseed" + "\n")

    fob.close()


    fob = open(Fold+"/steps.in","w+")

    fob.write("0.5d0\t,0.1d0\t,0.1d0\n")
    fob.write("step(1)\t,step(2)\t,step(3)\n")
    


   # fob.write("413000\t413000\t20\t200\t %s\t %s\n"%(str(G),str(random.randint(0,10E6)))) 
   # fob.write("runtime\tnsteps\tmappingsteps\ttrajectory\tG(eV)\tRandom\n")
   

def md(folder):
    try:
        os.mkdir(folder) 
    except:
        print "Folder exists:", folder
def cp(filename,folder):
    shutil.copy2(filename,folder)

nfold = 100

for i in range(nfold): 
    main = "mapp-" + str(i)
    md(main) 
    makeinp(main)
    cp('mat_1d.exe',main)
    cp('submit.sbatch',main)
    os.chdir(main)
    #sbatch("submit.sbatch")
    os.system("sbatch submit.sbatch")  
    os.chdir("../")
