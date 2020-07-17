import subprocess
import os

class SimulationEngine:

    def __init__(self, session):
        self._session = session

    def run(self):
        
        input_file = self._session._input_file
        output_file = self._session._output_file
        if self._session._command_type == "ev.x":
            command = ["ev.x"]
        else:
            command = [self._session._command_type, "-i", input_file, ">", output_file]
        try:
            proc = subprocess.run(" ".join(command), capture_output = True, shell = True)
        except:
            raise RuntimeError(f"An error occured when running the following command: {command}")