#!/usr/bin/python
import sys
import os
import subprocess
import time
import stat

def generate_file(in_filename, out_filename, params):
    out_ = open(out_filename,'w')
    with open(in_filename,'r') as in_:
        for line in in_:
            for p,v in params.items():
                line = line.replace(p, v)
            out_.write(line)
# {end def}

def check_queue(queue, n_max=4):
    n_running = 0
    queue2 = list(queue)
    queue = []
    for p in queue2:
        if not p.poll():
            n_running = n_running + 1
            queue.append(p)
    return len(queue) < n_max
# {end def}

def cleanup(string):
    return string.replace('.', '_')
# {end cleanup}

exec_dir = '/opt/development/polar_pfc/code/build'
out_dir = '/media/Daten/projects/polar_pfc'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

setup = 'vop'
run = 1

processes = [];
radius = 80
params = {'#{RADIUS}': str(radius)}

for r in [-0.98]:
    params['#{r}'] = str(r)

    for V0 in [0.3, 0.31, 0.32, 0.33, 0.34]:
        params['#{V0}'] = str(V0)
        params['#{POSTFIX}'] = cleanup('R' + str(radius) + '_r' + str(r) + '_v' + str(V0))
        print "postfix=",params['#{POSTFIX}']
        params['#{FILENAME}'] = 'solution_' + params['#{POSTFIX}']

        output = out_dir + '/results/' + setup + '/' + params['#{POSTFIX}']
        if not os.path.exists(output):
            os.makedirs(output)
        else:
            subprocess.call('rm -rf ' + output + '/*', shell=True)

        if not os.path.exists(output + '/data'):
            os.makedirs(output + '/data')

        params['#{OUT_DIR}'] = output
        params['#{INITFILE}'] = output + '/initfile'
        params['#{COMMAND}'] = exec_dir + '/polar_pfc ' + params['#{INITFILE}']

        generate_file('init/param.templ.json', params['#{INITFILE}'] + '.json', params)

        run_file = output + '/run.sh'
        generate_file('run.templ.sh', run_file, params)

        if os.path.exists(params['#{INITFILE}'] + '.json'):
            #command = 'sbatch ' + run_file
            command = params['#{COMMAND}']
            while not check_queue(processes):
                time.sleep(1)
            print command
            p = subprocess.Popen(command, shell=True)
            processes.append(p)

            print "Logfile will be written to ", params['#{OUT_DIR}'] + '/output.log'

        run = run+1
    # {end for}
# {end for}

