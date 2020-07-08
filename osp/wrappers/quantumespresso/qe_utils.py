from osp.core import QE
from osp.core.utils import simple_search
import re
import numpy as np  

class qeUtils():
    """Utilities for reading and writing .in and .out files
    """
    SystemSections = ('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS')


    def __init__(self, session, root):
        """__init__ function for using any of the following utils

        Args:
            session ([type]): the simulation CUDS object
        """
        self._session = session
        self._file_path_root = root
        self._simulation_type = "scf"

    def _create_input(self, file_name, sim, **kwargs):
        """Creates a .in file for QE. Enough information to run a simulation without **kwargs,
        but any real calculation will necessitate many.

        Args:
            file_path (str): file path of the .in file to be written 
            sim (CUDS object): CUDS simulation object
        """
        
        # Simulation parameters
        self.params = {
            "CONTROL": {
                "calculation": "'scf'",
                "pseudo_dir": "'.'",
                "tprnfor": ".true.",
                "tstress": ".true."
            },
            "SYSTEM": {
                "ibrav": 0,
                "ecutwfc": 100,
            },
                "ELECTRONS": {},
                "CELL": {},
        }
        # Updates simulation parameters based on keywords
        for key1, value1 in self.params.items():
            for key2, value2 in kwargs.items():
                if key1 == key2:
                    value1.update(value2)
        # Information about the system to be simulated
        self.sysinfo = {"ATOMIC_SPECIES":[], "ATOMIC_POSITIONS {crystal}":[],"K_POINTS {automatic}":[], "CELL_PARAMETERS {alat}":[]}
        
        def findo(oclass):
            return simple_search.find_cuds_objects_by_oclass(oclass = oclass, root = sim, rel = QE.HAS_PART)

        def _get_count(oclass):
            count = 0
            for item in findo(oclass):
                count += 1
            return count
        
        # Adds info from simulation CUDS to params (these are the only cases in which this is necessary)
        self.params["SYSTEM"]["nat"] = _get_count(QE.ATOM)
        self.params["SYSTEM"]["ntyp"] = _get_count(QE.ELEMENT)
        self.params["SYSTEM"]["celldm(1)"] = findo(QE.CELLDM1)[0].value
        
        self.atomlist = []

        # Adds info from simulation CUDS object to sysinfo dict
        for element in findo(QE.ELEMENT):
            self.sysinfo["ATOMIC_SPECIES"].append([element.name, element.get(oclass = QE.MASS)[0].value, element.get(oclass = QE.PSEUDOPOTENTIAL)[0].name])
        for atom in findo(QE.ATOM):
            self.atomlist.append(atom)
            print(atom.get(oclass = QE.POSITION)[0].vector)
            self.sysinfo["ATOMIC_POSITIONS {crystal}"].append([atom.get(oclass = QE.ELEMENT, rel = QE.IS_PART_OF)[0].name] + [atom.get(oclass = QE.POSITION)[0].vector[i] for i in range(3)])
        for point in findo(QE.K_POINTS):
            self.sysinfo["K_POINTS {automatic}"].append([int(point.vector[i]) for i in range(3)] + [0, 0, 0]) 
        for param in findo(QE.CELL_PARAMS)[0].get(rel = QE.HAS_PART):
            self.sysinfo["CELL_PARAMETERS {alat}"].append([i for i in param.vector])

        # Writes params to file
        with open(self._file_path_root+file_name, "w+") as f:
            for key1, value1 in self.params.items():
                f.write(f"&{key1} \n")
                for key2, value2 in value1.items():
                    f.write(f"  {key2} = {value2} \n")
                f.write("/\n")
            for key1, value1 in self.sysinfo.items():
                f.write(f" {key1} \n")
                for i in value1:
                    f.write(" ".join(str(v) for v in i) + "\n")

        # Gets the simulation type
        self._simulation_type = self.params["CONTROL"]["calculation"]

    def _read_output(self, file_name):
        """Reads the output file and returns energy

        Args:
            file_path (str): location of the output file

        Returns:
            energy: float that contains the energy of the system in Ry
        """
        with open(self._file_path_root+file_name, "r") as f:
            for line in f:
                if line.startswith("!"):
                    energy = float(line.split()[4])
        try: 
            return energy
        except:
            print("Error: output file did not contain total energy!")

    def update_cuds(self, file_name, sim):
        """Updates CUDS objects after the operation has been executed

        Args:
            simulation (Cuds): The simulation CUDS object
        """
        # TODO: Add stress tensor, also implement in ontology
        # TODO: Actually add/update these lol

        def findo(oclass):
            return simple_search.find_cuds_objects_by_oclass(oclass = oclass, root = sim, rel = QE.HAS_PART)

        with open(self._file_path_root+file_name, "r+") as file:
            pass
            
            # lines = file.readlines()
            # for i, line in enumerate(lines):
            #     if line.startswith("!"):
            #         total_energy = float(line.split()[4])
            #     if line.startswith("     atom "):
            #         atom = self.atomlist[int(line.split()[1])-1]
            #         if atom.get(oclass = QE.FORCE):
            #             atom.get(oclass = QE.FORCE)[0].vector = [float(line.split()[i]) for i in range(6, 9)]
            #         else:
            #             atom.add(QE.Force(vector = [float(line.split()[j]) for j in range(6, 9)], unit = ""))
            #         self.atomlist[int(line.split()[1])-1].get(oclass = QE.FORCE)[0].vector = [float(line.split()[i]) for i in range(6, 9)]
            #     if line.startswith("     Computing"):
            #         pressure = float(lines[i+2].split()[5])
            #         stresslines = [lines[i+j] for j in range(3, 6)]
            #         raw_stress_tensor = [float(i) for i in "".join(stresslines).split()]
            #         stress_tensor_kbar = np.array(raw_stress_tensor).reshape((3, 6))[:,3:6]

        print(self._simulation_type)

        def update_total_energy(line):
            if line.startswith("!"):
                total_energy = float(line.split()[4])
                cuds_entity = sim.get(oclass = QE.TotalEnergy)
                if cuds_entity:
                    cuds_entity[0].value = total_energy
                    cuds_entity[0].unit = "Ry"
                else:
                    sim.add(QE.TotalEnergy(value = total_energy, unit = "Ry"))

        def update_pressure(line):
            if line.startswith("     Computing"):
                pressure = float(lines[i+2].split()[5])
                cuds_entity = sim.get(oclass = QE.Pressure)
                if cuds_entity:
                    cuds_entity[0].value = pressure
                    cuds_entity[0].unit = "kbar"
                else:
                    sim.add(QE.Force(value = pressure, unit = "kbar"))

        def update_force(line):
            if line.startswith("     atom "):
                atom = self.atomlist[int(line.split()[1])-1]
                force = [float(line.split()[j]) for j in range(6, 9)]
                cuds_entity = atom.get(oclass = QE.Force)
                if cuds_entity:
                    cuds_entity[0].vector = force
                    cuds_entity[0].unit = "N"
                else:
                    atom.add(QE.Force(vector = force, unit = "N"))

        def update_stress_tensor(i, line):
            if line.startswith("     Computing"):
                stresslines = [lines[i+j] for j in range(3, 6)]
                raw_stress_tensor = [float(j) for j in "".join(stresslines).split()]
                stress_tensor_kbar = np.array(raw_stress_tensor).reshape((3, 6))[:,3:6]
                cuds_entity = sim.get(oclass = QE.StressTensor)
                if cuds_entity:
                    cuds_entity[0].tensor2 = stress_tensor_kbar
                    cuds_entity[0].unit = "kbar"
                else:
                    sim.add(QE.StressTensor(tensor2 = stress_tensor_kbar, unit = "kbar"))

        if self._simulation_type == "'scf'":
            with open(self._file_path_root+file_name, "r+") as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    update_total_energy(line)
                    update_pressure(line)
                    update_force(line)
                    update_stress_tensor(i, line)

        elif self._simulation_type == "'relax'":
            with open(self._file_path_root+file_name, "r+") as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    update_total_energy(line)
                    update_pressure(line)
                    update_force(line)

            # Update atomic forces, pressure, total energy, stress tensor

        # elif self._simulation_type == 'relax':
        #     # Update atomic forces, pressure, total energy, stress tensor, final atomic coordinates
        
        # elif self._simulation_type == 'vc-relax':
        #     # Update atomic forces, pressure, total energy, stress tensor, final atomic coordinates, final cell parameters
