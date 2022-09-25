import itertools
import os
import sys
import subprocess
import time

#IMPORT MODULES FROM OTHER DIR

Geom_1 = ["D92"]
Shift = ["1", "2", "3", "4", "5", "6", "7"]

if not os.path.exists("tmpSub1/log"):
    os.makedirs("tmpSub1/log")
condorLogDir = "log"
tarFile = "tmpSub1/generator.tar.gz"
if os.path.exists(tarFile):
    os.system("rm %s"%tarFile)
os.system("tar -zcvf %s ../../relval --exclude condor"%tarFile)
os.system("cp rungen.sh tmpSub1/")
common_command = \
'Universe   = vanilla\n\
should_transfer_files = YES\n\
when_to_transfer_output = ON_EXIT\n\
Transfer_Input_Files = generator.tar.gz, rungen.sh\n\
x509userproxy = $ENV(X509_USER_PROXY)\n\
use_x509userproxy = true\n\
RequestCpus = 8\n\
+BenchmarkJob = True\n\
#+JobFlavour = "testmatch"\n\
#+MaxRuntime = 259200\n\
+MaxRuntime = 7200\n\
Output = %s/log_$(cluster)_$(process).stdout\n\
Error  = %s/log_$(cluster)_$(process).stderr\n\
Log    = %s/log_$(cluster)_$(process).condor\n\n'%(condorLogDir, condorLogDir, condorLogDir)

#----------------------------------------
#Create jdl files
#----------------------------------------
geom_par = 1
sampleList = eval("Geom_%i"%(geom_par))
shiftList = eval("Shift")
jdlName = 'submitJobs_%s.jdl'%(geom_par)
jdlFile = open('tmpSub1/%s'%jdlName,'w')
jdlFile.write('Executable =  rungen.sh \n')
jdlFile.write(common_command)
jdlFile.write("X=$(step)\n")
for sample in sampleList:
    for shift in shiftList:
        condorOutDir1="/eos/cms/store/group/dpg_hgcal/comm_hgcal/geomval/cassette_shift_mult"
        #condorOutDir1="/eos/user/i/idas/SimOut/geomval/muontomo_newtrig"
        os.system("eos root://eosuser.cern.ch mkdir -p %s/%s/%s"%(condorOutDir1, sample, shift))
        condorOutDir="/cms/store/user/idas/SimOut/geomval/cassette_shift_mult"
        os.system("xrdfs root://se01.indiacms.res.in/ mkdir -p %s/%s/%s"%(condorOutDir, sample, shift))
        run_command =  'Arguments  = %s %s $INT(X) \nQueue 10\n\n' %(sample, shift)
        jdlFile.write(run_command)
        #print "condor_submit jdl/%s"%jdlFile
jdlFile.close() 

