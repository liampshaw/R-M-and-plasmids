# Removes directories - but checks if the request is dangerous for safety
import os
import sys
import shutil
import argparse

def get_options():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dir', nargs='+', help='directory to remove', required=False, default='')
    return parser.parse_args()

args = get_options()
dir = args.dir[0] 
if len(dir)>4: # don't trust short arguments
    if not '*' in dir: # don't trust wildcards
        if os.path.exists(dir):
            shutil.rmtree(dir)


