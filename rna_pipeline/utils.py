# utils.py
# Group Members: Moe Sithu Maung Maung Lay, Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh

'''
utils.py

Purpose: Provides the pipeline with useful tools, helpers and anything important that may be used anywhere else.
'''

#Imports path for filesystem paths, subprocess for external programs, 
#Also imports sys and logging for Python's logging system, yaml to load yml config files
from pathlib import Path
import subprocess
import sys
import logging
import yaml

def getLogger (name):
    '''
    Function: getLogger
    Purpose: Gives every module a consistent logging style, and makes them go to the same place
    - Prints labeled messages to the console, and avoids double-adding handlers if called multiple times.
    Inputs: name (str)
    Outputs: logging.Logger object
    '''
    
    #Creates a logger with a given name
    logger = logging.getLogger(name)
    
    #Prevents adding duplicate handlers if it gets called multiple times
    #Applies a log format and returns the configured logger
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(levelname)s] %(name)s: %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger

def runCmd (cmd, logFile=None):
    '''
    Function: runCmd
    Purpose: Allows us to run external tools in the pipeline like STAR, BLAST, FastQC
    - Everything goes through cmd so behavior is consistent and errors are caught in the same way
    Inputs: cmd (list of str), logFile (Path or None)
    Outputs: None
    '''

    #Prints the message for the user to see what they are running and what's happening
    print(f"Running command: {' '.join(cmd)}")
    
    #Runs the command with or without logging to file
    #If provided, it opens the file for writing, stdout and stderr go in there
    #If not provided, it just runs and it goes into the terminal
    if logFile:
        with logFile.open("w") as lf:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    else:
        proc = subprocess.run(cmd)

    #If commandline fails it shows an error   
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


#Helper that creates directories and checks for parent directories as well
#If they already exist then it does not do anything
def ensureDir (path):
    '''
    Create directory if it does not exist.
    Inputs: path (Path)
    Outputs: None
    '''
    path.mkdir(parents=True, exist_ok=True)

#Reads and loads a YAML reference configuration file and returns it as a Python dictionary
def loadReferenceConfig (path):
    '''
    Load YAML reference configuration file.
    Inputs: path (Path or None)
    Outputs: dict
    '''
    
    #If no config given, returns an empty dictionary
    if path is None:
        return {}
    
    #Handle string path input safely
    if isinstance(path, str):
        path = Path(path)

    #Makes it into python dictionary  
    with path.open() as f:
        return yaml.safe_load(f)

#Function to create and return a named subdirectory under a given base directory
def subDir (base, name):
    '''
    Create and return a subdirectory of base.
    Inputs: base (Path), name (str)
    Outputs: Path
    '''
    p = base / name
    ensureDir(p)
    return p