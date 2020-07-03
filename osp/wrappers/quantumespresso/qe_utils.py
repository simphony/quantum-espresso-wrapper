from osp.core import QE
from osp.core.utils import simple_search

class qeUtils():
    """Utilities for reading and writing .in and .out files
    """
    SystemSections = ('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS')


    def __init__(self, session):
        """__init__ function for using any of the following utils

        Args:
            session ([type]): the simulation CUDS object
        """
        self._session = session

    def _create_input(self, file_path, sim, **kwargs):
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

        # adds info from simulation CUDS object to sysinfo dict
        for element in findo(QE.ELEMENT):
            self.sysinfo["ATOMIC_SPECIES"].append([element.name, element.get(oclass = QE.MASS)[0].value, element.get(oclass = QE.PSEUDOPOTENTIAL)[0].name])
        for atom in findo(QE.ATOM):
            print(atom.get(oclass = QE.POSITION)[0].vector)
            self.sysinfo["ATOMIC_POSITIONS {crystal}"].append([atom.get(oclass = QE.ELEMENT, rel = QE.IS_PART_OF)[0].name] + [atom.get(oclass = QE.POSITION)[0].vector[i] for i in range(3)])
        for point in findo(QE.K_POINTS):
            self.sysinfo["K_POINTS {automatic}"].append([int(point.vector[i]) for i in range(3)] + [0, 0, 0]) 
        for param in findo(QE.CELL_PARAMS)[0].get(rel = QE.HAS_PART):
            self.sysinfo["CELL_PARAMETERS {alat}"].append([i for i in param.vector])

        # Writes params to file
        with open(file_path, "w+") as f:
            for key1, value1 in self.params.items():
                f.write(f"&{key1} \n")
                for key2, value2 in value1.items():
                    f.write(f"  {key2} = {value2} \n")
                f.write("/\n")
            for key1, value1 in self.sysinfo.items():
                f.write(f" {key1} \n")
                for i in value1:
                    f.write(" ".join(str(v) for v in i) + "\n")

    def _read_input(self, file_path):
        # TODO: finish this and decide on a concrete way for this to work
        """Reads input to create CUDS structure automatically

        Args:
            file_path (str): location of the .in file

        Returns:
            TBD: TBD
        """
        output = {}
        with open(file_path, "r") as f:
            text = _read_clean(f)
            for line in f:
                if line.startswith(self.SystemSections):
                    section = line
                    output[section] = []
                else:
                    if line != '/':
                        output[section].append(line)
        return output

    def _read_clean(self, file):
        content = []
        for line in file.readlines():
            line = line.partition('#')[0].strip()
            if line:
                content.append(line)
        return content

    def _read_output(self, file_path):
        """Reads the output file and returns energy

        Args:
            file_path (str): location of the output file

        Returns:
            energy: float that contains the energy of the system in Ry
        """
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("!"):
                    energy = float(line.split()[4])
        try: 
            return energy
        except:
            print("Error: output file did not contain total energy!")