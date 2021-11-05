
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import time
import yaml
import argparse
import subprocess
import numpy as np

# color
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\33[93m'

def setup():
  for phase in phases:
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      with open(fname + ".key", 'w') as fw:
        for line in phase_key[phase]:
          if 'parameters' in line.lower():
            line = f'parameters     ../{prm}\n'
          fw.write(line)
        fw.write('\n')
        fw.write(f'ligand -1 {natom}\n')
        fw.write(f'ele-lambda {elb}\n')
        fw.write(f'vdw-lambda {vlb}\n')
      linkxyz = f"ln -s {homedir}/{phase_xyz[phase]} {homedir}/{phase}/{fname}.xyz"
      movekey = f"mv {fname}.key ./{phase}"
      subprocess.run(linkxyz, shell=True)
      subprocess.run(movekey, shell=True)
  print(GREEN + "BAR simulation files generated!" + ENDC)
  return

def dynamic():
  phase_dynamic = {'liquid':'/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2/dynamic9.sh', \
                   'gas': '/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8/dynamic.x'}
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      xyzfile = fname + ".xyz"
      keyfile = xyzfile.replace("xyz", "key")
      logfile = xyzfile.replace("xyz", "log")
      if (os.path.isfile(keyfile)) and (not os.path.isfile(logfile)):
        if phase == 'liquid':
          if liquidensemble == "NPT":
            dynamiccmd = f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 4 {liquidtemperature} {liquidpressure} > {logfile}"
          elif liquidensemble == "NVT":
            dynamiccmd = f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 2 {liquidtemperature} > {logfile}"
          else:
            sys.exit(RED + "Error: only NPT or NVT ensemble is supported for Free energy simulations" + ENDC)
          if nodes == []:
            subprocess.run(dynamiccmd, shell=True)
          else:
            subprocess.run(f"nohup python {submitexe} -c '{dynamiccmd}' -n {nodes[i]} 2>err &", shell=True)
        if phase == 'gas':
          dynamiccmd = f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {gastotalstep} {gastimestep} {gaswriteout} 2 {gastemperature} > {logfile}"
          if nodes == []:
            subprocess.run(dynamiccmd, shell=True)
          else:
            subprocess.run(f"nohup python {submitexe} -c '{dynamiccmd}' -n {nodes[i]} 2>err &", shell=True)
  return

def bar():
  proceed = True
  phase_simtime = {'liquid':1000.0*liquidtotaltime, 'gas':1000.0*gastotaltime}
  print(YELLOW + "\n Checking the completeness of the MD trajectories" + ENDC)
  for phase in phases:
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      logfile = fname + ".log"
      if os.path.isfile(os.path.join(homedir, phase, logfile)):
        lines = open(os.path.join(homedir, phase, logfile)).readlines()
        for line in lines[-20:]:
          if "Current Time" in line:
            simtime = float(line.split()[2])
            if (simtime == phase_simtime[phase]):
              print(GREEN + "  [" + fname + f"]: finished {simtime} out of {phase_simtime[phase]} ps!" + ENDC)
            else:
              print(RED + "  [" + fname + f"]: finished {simtime} out of {phase_simtime[phase]} ps!" + ENDC)
              proceed = False
  if proceed:
    barexe = '/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8/bar.x'
    for phase in phases:
      for i in range(0,len(orderparams)-1,1):
        elb0, vlb0 = orderparams[i]
        elb1, vlb1 = orderparams[i+1]
        fname0 = "%s-e%s-v%s"%(phase, "%03d"%int(elb0*100), "%03d"%int(vlb0*100))
        fname1 = "%s-e%s-v%s"%(phase, "%03d"%int(elb1*100), "%03d"%int(vlb1*100))
        arcfile0 = fname0 + ".arc"
        arcfile1 = fname1 + ".arc"
        if phase == 'liquid':
          outfile = fname0 + ".out"
          barfile = fname0 + ".bar"
          enefile = fname0 + ".ene"
          barcmd1 = f"{barexe} 1 {arcfile0} {liquidtemperature} {arcfile1} {liquidtemperature} N > {outfile} && "
          totalsnapshot = int(phase_simtime[phase]/liquidwriteout)
          startsnapshot = int(totalsnapshot/5.0) + 1
          barcmd2 = f"{barexe} 2 {barfile} {startsnapshot} {totalsnapshot} 1 {startsnapshot} {totalsnapshot} 1 > {enefile} "
          barstr = barcmd1 + barcmd2
          if nodes == []:
            subprocess.run(cmdstr, shell=True)
          else:
            phasedir = os.path.join(homedir, phase)
            os.chdir(phasedir)
            cmdstr = f"nohup python {submitexe} -c '{barstr}' -n {nodes[i]} 2>err &"
            subprocess.run(cmdstr, shell=True)
            print(GREEN + f"submitted {barstr} on {nodes[i]}" + ENDC)
        if phase == 'gas':
          outfile = fname1 + ".out"
          barfile = fname1 + ".bar"
          enefile = fname1 + ".ene"
          barcmd1 = f"{barexe} 1 {arcfile1} {gastemperature} {arcfile0} {gastemperature} N > {outfile} && "
          totalsnapshot = int(phase_simtime[phase]/liquidwriteout)
          startsnapshot = int(totalsnapshot/5.0) + 1
          barcmd2 = f"{barexe} 2 {barfile} {startsnapshot} {totalsnapshot} 1 {startsnapshot} {totalsnapshot} 1 > {enefile} "
          barstr = barcmd1 + barcmd2
          if nodes == []:
            subprocess.run(barstr, shell=True)
          else:
            phasedir = os.path.join(homedir, phase)
            os.chdir(phasedir)
            cmdstr = f"nohup python {submitexe} -c '{barstr}' -n {nodes[i]} 2>err &"
            subprocess.run(cmdstr, shell=True)
            print(GREEN + f"submitted {barstr} on {nodes[i]}" + ENDC)
  return

def result():
  proceed = True
  gasenes = []
  liquidenes = []
  gasperturbsteps = []
  liquidperturbsteps = []
  for phase in phases:
    for i in range(0,len(orderparams)-1,1):
      elb0, vlb0 = orderparams[i]
      elb1, vlb1 = orderparams[i+1]
      fname0 = "%s-e%s-v%s"%(phase, "%03d"%int(elb0*100), "%03d"%int(vlb0*100))
      fname1 = "%s-e%s-v%s"%(phase, "%03d"%int(elb1*100), "%03d"%int(vlb1*100))
      if phase == 'gas':
        enefile = fname1 + ".ene"
        gasenes.append(enefile)
        gasperturbsteps.append([fname1, fname0])
      if phase == 'liquid':
        enefile = fname0 + ".ene"
        liquidenes.append(enefile)
        liquidperturbsteps.append([fname0, fname1])
      if not os.path.isfile(os.path.join(homedir, phase, enefile)):
        print(RED + fname + f": free energy file (.ene) not found!" + ENDC)
        proceed = False
  if proceed:
    FEgas = []
    Errgas = []
    FEliquid = []
    Errliquid = []
    for gasene in gasenes:
      find = False
      for line in open(os.path.join('gas', gasene)).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          FEgas.append(float(line.split()[-4]))
          Errgas.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          FEgas.append(float(line.split()[-4]))
          Errgas.append(float(line.split()[-2]))
      if not find:
        sys.exit(RED + f"Could not find free energy for {gasene}" + ENDC)
    
    for liquidene in liquidenes:
      find = False
      for line in open(os.path.join('liquid', liquidene)).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          FEliquid.append(float(line.split()[-4]))
          Errliquid.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          FEliquid.append(float(line.split()[-4]))
          Errliquid.append(float(line.split()[-2]))
      if not find:
        sys.exit(RED + f"Could not find free energy for {liquidene}" + ENDC)
    
    FEall = FEgas + FEliquid
    Errall = Errgas + Errliquid
    totFE = np.array(FEall).sum()
    totErr = np.sqrt(np.square(np.array(Errall)).sum())
    fo = open("result.txt", "w")
    print(GREEN + "%20s%20s%30s%20s"%("StateA", "StateB", "FreeEnergy(kcal/mol)", "Error(kcal/mol)") + ENDC)
    fo.write("%20s%20s%30s%20s\n"%("StateA", "StateB", "FreeEnergy(kcal/mol)", "Error(kcal/mol)"))
    for i in range(len(FEgas)-1,-1,-1):
      state0, state1 = gasperturbsteps[i]
      print("%20s%20s%25.4f%20.4f"%(state0, state1, FEgas[i], Errgas[i]))
      fo.write("%20s%20s%25.4f%20.4f\n"%(state0, state1, FEgas[i], Errgas[i]))
    for i in range(len(FEliquid)):
      state0, state1 = liquidperturbsteps[i]
      print("%20s%20s%25.4f%20.4f"%(state0, state1, FEliquid[i], Errliquid[i]))
      fo.write("%20s%20s%25.4f%20.4f\n"%(state0, state1, FEliquid[i], Errliquid[i]))
    print(GREEN + "%40s%25.4f%20.4f"%("SUMMARY OF THE TOTAL FREE ENERGY", totFE, totErr) + ENDC)
    fo.write("%40s%25.4f%20.4f\n"%("SUMMARY OF THE TOTAL FREE ENERGY", totFE, totErr))
    fo.close()
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument( 'act', help = "Actions to take.", choices = ['setup', 'dynamic', 'bar', 'result'])  
  args = vars(parser.parse_args())
  
  global rootdir
  rootdir = os.path.join(os.path.split(__file__)[0])

  if not os.path.isfile("settings.yaml"):
    sys.exit(RED + "Please provide 'settings.yaml' file; " + ENDC + GREEN + f"An example is here for you {os.path.join(rootdir, 'dat', 'settings.yaml')}" + ENDC)
  else:
    with open('settings.yaml') as f:
      FEsimsettings = yaml.load(f, Loader=yaml.FullLoader)
  
  global prm, run
  lig = FEsimsettings['gas_xyz'] 
  box = FEsimsettings['box_xyz'] 
  prm = FEsimsettings["tinker_prm"]
  hfe = FEsimsettings["hydration"] 
  nodelist = FEsimsettings["node_list"]
 
  global homedir
  homedir = os.getcwd()

  global orderparams, liquidkeylines, gaskeylines
  orderparams = []
  orderprmfile = os.path.join(rootdir, "dat", "orderparams")
  liquidkeylines = open(os.path.join(rootdir, "dat", "liquid.key")).readlines()
  gaskeylines = open(os.path.join(rootdir, "dat", "gas.key")).readlines()
  
  for line in open(orderprmfile).readlines():
    if "#" not in line:
      d = line.split()
      orderparams.append([float(d[0]), float(d[1])])
  
  global phases
  if hfe:
    phases = ['liquid', 'gas']
  else:
    phases = ['liquid']
  for phase in phases:
    if not os.path.isdir(f"{phase}"):
      os.system(f"mkdir {phase}")
  
  global submitexe
  submitexe = os.path.join(rootdir, "utils", "submit.py")
  
  global nodes
  nodes = FEsimsettings["node_list"]
  if nodes == None:
    nodes = []

  global natom, phase_xyz, phase_key
  natom = int(open(lig).readlines()[0].split()[0])
  phase_xyz = {'liquid':box, 'gas':lig}
  phase_key = {'liquid':liquidkeylines, 'gas':gaskeylines}

  global liquidtotaltime, liquidtimestep, liquidwriteout, liquidtemperature, liquidpressure, liquidensemble
  global liquidtotalstep
  liquidtotaltime = FEsimsettings["liquid_md_total_time"] 
  liquidtimestep = FEsimsettings["liquid_md_time_step"] 
  liquidwriteout = FEsimsettings["liquid_md_write_freq"] 
  liquidtemperature = FEsimsettings["liquid_md_temperature"] 
  liquidpressure = FEsimsettings["liquid_md_pressure"]
  liquidensemble = FEsimsettings["liquid_md_ensemble"].upper()
  liquidtotalstep = int((1000000.0*liquidtotaltime)/liquidtimestep)

  global gastotaltime, gastimestep, gaswriteout, gastemperature, gastotalstep
  gastotaltime = FEsimsettings["gas_md_total_time"] 
  gastimestep = FEsimsettings["gas_md_time_step"] 
  gastemperature = FEsimsettings["gas_md_temperature"] 
  gaswriteout = FEsimsettings["gas_md_write_freq"]
  gastotalstep = int((1000000.0*gastotaltime)/gastimestep)

  actions = {'setup':setup, 'dynamic':dynamic, 'bar':bar, 'result':result}
  actions[args['act']]()
  return

if __name__ == "__main__":
  main()
