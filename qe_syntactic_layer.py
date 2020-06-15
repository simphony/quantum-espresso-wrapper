import subprocess

def run_pw(inputfilename):
    subprocess.run(['pw.x', '-input', inputfilename, '> test.out'])

