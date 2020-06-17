import os
import logging
import subprocess

class SimulationEngine:

    def __init__(self, session):
        self._session = session

    def run(self):
        """Run the energy calculation
        """
        self._execute_command(
            arguments = [],
            input_files={
                "-input": "inputfile.in"
            },
            output_files = {
                ">": "basicSiCrystal.out"
            })

    def _execute_command(self, arguments, input_files, output_files):
        """Executes a command in console

        Args:
            arguments (str): any additional arguments that may be needed
            input_files (str): name of the file that contains input data
            output_files (str): name of the file that will contain output data

        Raises:
            RuntimeError: An error occured when executing the command.
        """
        command = ["pw.x"] + arguments
        command += [x for item in input_files.items() for x in item]
        command += [x for item in output_files.items() for x in item]
        print(" ".join(command))
        try:
            proc = subprocess.run(" ".join(command), capture_output = True, shell = True)
        except subprocess.CalledProcessError:
            raise RuntimeError("An error occured when running %s" % command)