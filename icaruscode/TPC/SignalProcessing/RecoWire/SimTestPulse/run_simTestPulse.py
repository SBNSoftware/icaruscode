import sys,os,commands,subprocess,time

if len(commands.getoutput('which lar').split())>1:
    sys.stderr.write('No "lar" command found...\n')
    sys.stderr.flush()
    sys.exit(1)

WORK_DIR='/tmp/%s_fakenews' % os.environ['USER']
OUT_FILE ='simTestPulse_full.root'

if os.path.isdir(WORK_DIR):
    sys.stderr.write('Found already-exisitng work dir "%s" ... remove and try re-run!\n' % WORK_DIR)
    sys.stderr.flush()
    sys.exit(1)
if os.path.isfile(OUT_FILE):
    sys.stderr.write('Found already-exisitng out file "%s" ... remove and try re-run!\n' % OUT_FILE)
    sys.stderr.flush()
    sys.exit(1)    

os.mkdir(WORK_DIR)

proc  = subprocess.Popen('cd %s;lar -c simtestpulse_driver.fcl' % WORK_DIR,shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
state = None
time_slept = 0
while state is None:
    state = proc.poll()
    sys.stdout.write('Running lar... (%-2d [s])\r' % time_slept)
    sys.stdout.flush()
    time.sleep(1)
    time_slept += 1
print

fout=open('log.txt','w')
out,err = proc.communicate()
fout.write(out)
fout.close()
print 'Log file created: log.txt'

if not state == 0:
    sys.stderr.write('Failed running the code (state=%s)! Check log.txt...\n' % str(state))
    sys.stderr.flush()
else:
    os.system('hadd %s %s/simTestPulse.root %s/simTestPulseAna.root' % (OUT_FILE,WORK_DIR,WORK_DIR))
    print
    print 'Output:',OUT_FILE
    print
os.system('rm -r %s' % WORK_DIR)
