import sys
import argparse

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+', help='command to put in qsub script')
    parser.add_argument('--name', help='name of job', required=False, default='job')
    parser.add_argument('--array', help='number of arrays', required=False, default='none')
    return parser.parse_args()

args = get_options()
script_command=str(args.input[0])
name = str(args.name)
array = str(args.array)
print("#!/bin/bash")

print("#$ -cwd -V")
print("#$ -N %s -j y" % name)
print("#$ -q short.qc")
if array!='none':
    print("#$ -t %s" % array)

print("echo \"****************************************************\"")
print("echo \"SGE job ID: \"$JOBID")
print("echo \"SGE task ID: \"$SGE_TASK_ID")
print("echo \"Run on host: \"`hostname`")
print("echo \"Operating system: \"`uname -s`")
print("echo \"Username: \"`whoami`")
print("echo \"Started at: \"`date`")
print("echo \"****************************************************\"")

print("# CONDA command")
print("module load Anaconda3/2020.11")
print("eval \"$(conda shell.bash hook)\"")


print("# COMMANDS GO HERE")
print(script_command)


print("echo \"****************************************************\"")
print("echo \"Finished at: \"`date`")
print("echo \"****************************************************\"")
