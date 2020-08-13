import subprocess
import pexpect
import osp.wrappers.quantumespresso.qe_utils

class SimulationEngine:
    def __init__(self, session):
        self._session = session

    def run(self):
        
        input_file = self._session._qe_utils._file_path_root + self._session._input_file
        output_file = self._session._qe_utils._file_path_root + self._session._output_file

        # Using pexpect, interacts with the ev.x command using certain variables
        # Couldn't find any other way to do this, if someone can do it better, please let me know
        if self._session._qe_utils.__class__.__base__ == osp.wrappers.quantumespresso.qe_utils.cliUtils:
            child = pexpect.spawn(self._session._command_type)
            for i, j in self._session._qe_utils.params.items():
                child.expect(i)
                child.sendline(j)
            child.wait()

        # if self._session._command_type == "ev.x":
        #     child = pexpect.spawn('ev.x')
        #     child.expect('au')
        #     child.sendline('au')
        #     child.expect('fcc')
        #     child.sendline('noncubic')
        #     child.expect('1')
        #     child.sendline(self._session._calculation_type)
        #     child.expect('Input')
        #     child.sendline(input_file)
        #     child.expect('Output')
        #     child.sendline(output_file)
        #     child.wait()

        # Runs the command in the usual way
        else:
            command = [self._session._command_type, "-i", input_file, ">", output_file]
            try:
                proc = subprocess.run(" ".join(command), capture_output = True, shell = True)
                print(" ".join(command))
            except:
                raise RuntimeError(f"An error occured when running the following command: {command}")