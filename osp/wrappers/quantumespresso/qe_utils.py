from osp.core.namespaces import QE
from osp.core.utils import simple_search
import numpy as np  

class qeUtils():
    """Utilities for reading and writing .in and .out files
    """


    def __init__(self, session, root):
        """__init__ function for using any of the following utils

        Args:
            session (cuds object): the simulation CUDS object
        """
        self._session = session
        self._file_path_root = root
        self.params = {}

    def _create_input(self, sim, **kwargs):
        """Creates input file(s) necessary to perform the calculations

        Args:
            sim (QE.Simulation or list of QE.Simulations): the simulation on which to perform the calculation.
            For calculations that require multiple simulations and aggregate the data (such as ev.x), please provide a list of strings.

        """
        # Update params based on kwargs
        for key1, value1 in self.params.items():
            for key2, value2 in kwargs.items():
                if key1 == key2:
                    value1.update(value2)

        # Writes to file based on params and sys
        print(self._file_path_root + self._session._input_file)
        with open(self._file_path_root + self._session._input_file, "w+") as f:
            for key1, value1 in self.params.items():
                f.write(f"&{key1} \n")
                for key2, value2 in value1.items():
                    f.write(f"  {key2} = {value2} \n")
                f.write("/\n")
            if self._session._command_type == "pw.x" or self._session._command_type == "ph.x": # TODO: find a way to put this in the pwUtils class
                for key1, value1 in self.sysinfo.items():
                    f.write(f" {key1} ")
                    for i in value1:
                        f.write(" ".join(str(v) for v in i) + "\n")

    def _update_cuds(self, sim):
        """Based off of the structure

        Args:
            sim (QE.Simulation or list of QE.Simulations): the simulation for which cuds should be updated.
            For calculations that require multiple simulations and aggregate the data (such as ev.x), please provide a list of strings.
        """
        pass

class pwUtils(qeUtils):
    def _create_input(self, sim, **kwargs):
        # Simulation parameters
        self.params = {
                "CONTROL": {
                    "calculation": f"'{self._session._calculation_type}'",
                    "pseudo_dir": "'.'",
                    "tprnfor": ".true.",
                    "tstress": ".true.",
                    "prefix": f"'{self._session._prefix}'",
                },
                "SYSTEM": {
                    "ibrav": 0,
                    "ecutwfc": 100,
                },
                    "ELECTRONS": {},
                    "CELL": {},
                    "IONS": {}
            }

        # Information about the system to be simulated
        self.sysinfo = {"ATOMIC_SPECIES":[[""]], "ATOMIC_POSITIONS":[["{crystal}"]],"K_POINTS":[["{automatic}"]], "CELL_PARAMETERS":[["{alat}"]]}

        # Defining a couple useful functions
        def _get_count(oclass):
            count = 0
            for i in simple_search.find_cuds_objects_by_oclass(oclass = oclass, root = sim, rel = QE.HAS_PART):
                count +=1
            return count

        def findo(oclass, depth):
            return simple_search.find_cuds_objects_by_oclass(oclass = oclass, root = sim, rel = QE.HAS_PART)

        # Add some sysinfo based on cuds
        self.params["SYSTEM"]["nat"] = _get_count(oclass = QE.Atom)
        self.params["SYSTEM"]["ntyp"] = _get_count(QE.Element)
        self.params["SYSTEM"]["celldm(1)"] = findo(QE.Celldm1, 2)[0].value

        # Storing atoms so that the same order can be used to update cuds later on
        self.atomlist = []

        # Adds a bunch of stuff to sysinfo
        for element in findo(QE.Element, 1):
            self.sysinfo["ATOMIC_SPECIES"].append([element.name, element.get(oclass = QE.Mass)[0].value, element.get(oclass = QE.PSEUDOPOTENTIAL)[0].name])
        
        for atom in findo(QE.Atom, 2):
            self.atomlist.append(atom)
            self.sysinfo["ATOMIC_POSITIONS"].append([atom.get(oclass = QE.Element, rel = QE.IS_PART_OF)[0].name] + [i for i in atom.get(oclass = QE.Position)[0].vector])

        for point in findo(QE.K_POINTS, 1):
            self.sysinfo["K_POINTS"].append([int(i) for i in point.vector] + [0, 0, 0])

        for param in findo(QE.CellParams, 2)[0].get(rel = QE.HAS_PART):
            self.sysinfo["CELL_PARAMETERS"].append([i for i in param.vector])
        
        # Inherits method
        super()._create_input(sim, **kwargs)

    def _update_cuds(self, sim):
        # Adds the output file in case anyone wants to take a look at it
        sim.add(QE.PwOut(path = self._file_path_root + self._session._output_file))
            
        # A variety of functions to update particular aspects of a cuds simulation
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

        def update_atomic_positions(i, line):
            if line.startswith("Begin"):
                positionslines = [lines[i+j] for j in range(3, 3+len(self.atomlist))]
                for j, line in enumerate(positionslines):
                    atom = self.atomlist[j]
                    position = [float(line.split()[k]) for k in range(1, 4)]
                    cuds_entity = atom.get(oclass = QE.Position)
                    cuds_entity[0].vector = position
                    cuds_entity[0].unit = "kbar"

        def update_celldm1(line):
            if line.startswith("CELL_PARAMETERS"):
                celldm1 = float(line.split()[2][:-1])
                cuds_entity = sim.get(oclass = QE.Cell)[0].get(oclass = QE.Celldm1)
                cuds_entity[0].value = celldm1
                cuds_entity[0].unit = "au"

        def update_cell_params(i, line):
            if line.startswith("CELL_PARAMETERS"):
                paramlines = [lines[i+j] for j in range(1, 4)]
                cuds_entity = sim.get(oclass = QE.Cell)[0].get(oclass = QE.Cell)[0].get(oclass = QE.CellParams)[0]
                cuds_entity.get(oclass = QE.CellParameterX)[0].vector = [float(k) for k in paramlines[0].split()]
                cuds_entity.get(oclass = QE.CellParameterY)[0].vector = [float(k) for k in paramlines[1].split()]
                cuds_entity.get(oclass = QE.CellParameterZ)[0].vector = [float(k) for k in paramlines[2].split()]

        def update_volume(line):
            if line.startswith("     unit-cell volume"):
                volume = float(line.split()[3])
                cuds_entity = sim.get(oclass = QE.Cell)[0].get(oclass = QE.Volume)
                if cuds_entity:
                    cuds_entity[0].value = volume
                    cuds_entity[0].unit = "au^3"
                else:
                    sim.get(oclass = QE.Cell)[0].add(QE.Volume(value = volume, unit = "au^3"))

        # How the cuds simulation should be updated depending on what calculation type
        # Test
        if self._session._calculation_type == "scf":
            with open(self._file_path_root + self._session._output_file, "r+") as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    update_total_energy(line)
                    update_pressure(line)
                    update_force(line)
                    update_stress_tensor(i, line)
                    update_volume(line)

        if self._session._calculation_type == "relax":
            with open(self._file_path_root + self._session._output_file, "r+") as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    update_total_energy(line)
                    update_pressure(line)
                    update_force(line)
                    update_stress_tensor(i, line)
                    update_atomic_positions(i, line)
                    update_volume(line)

        
        if self._session._calculation_type == "vc-relax":
            with open(self._file_path_root + self._session._output_file, "r+") as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    update_total_energy(lines)
                    update_pressure(lines)
                    update_force(line)
                    update_stress_tensor(i, line)
                    update_atomic_positions(i, line)
                    update_celldm1(i, line)
                    update_volume(line)
        super()._update_cuds(sim)

class bandsUtils(qeUtils):
    def _create_input(self, sim, **kwargs):
        self.params = {
                "BANDS": {
                    "prefix": f"'{self._session._prefix}'",
                    "outdir": "'.'",
                    "filband": f"'{self._session._prefix}" + ".bands.dat'"
                }
            }
        self._session._calculation_type = ""
        self.sysinfo = []
        super()._create_input(sim, **kwargs)

    def _update_cuds(self, sim):
        sim.add(QE.BandsDat(path = self._file_path_root + self._session._prefix + ".bands.dat"))
        super()._update_cuds(sim)

class dosUtils(qeUtils):
    def _create_input(self, sim, **kwargs):
        self.params = {
                "DOS": {
                    "outdir": "'.'",
                    "prefix": f"'{self._session._prefix}'",
                    "DeltaE": 0.05,
                    "fildos": f"'{self._session._prefix}" + ".dos.dat'"
                }
            }
        super()._create_input(sim, **kwargs)

    def _update_cuds(self, sim):
        sim.add(QE.DosDat(path = self._file_path_root + self._session._prefix + ".dos.dat"))
        super()._update_cuds(sim)

class ppUtils(qeUtils):
    def _create_input(self, sim, **kwargs):
        self.params = {
                "INPUTPP": {
                    "prefix": f"'{self._session._prefix}'",
                    "outdir":  "'.'",
                    "filplot": f"'{self._session._prefix}.pp{self._session._calculation_type}.txt'",
                    # Note that plot_num is strictly int, reference to the significance of each values can be found here: https://www.quantum-espresso.org/Doc/INPUT_PP.html
                    # We use calculation type because it is already in use
                    "plot_num": self._session._calculation_type
                },
                "PLOT": {
                    "nfile": 1,
                    "iflag": 3,
                    "output_format": 3,
                    "fileout": f"'{self._session._prefix}{self._session._calculation_type}.pp.xsf'",
                    "nx": 101,
                    "ny": 101,
                    "nz": 101
                    # TODO: add support for manual vectors here
                    # TODO: add support for variable output formats
                }

            }
        super()._create_input(sim, **kwargs)

    def _update_cuds(self, sim):
        sim.add(QE.XSF(path = self._file_path_root + self._session._prefix + ".pp.xsf"))
        super()._update_cuds(sim)

class evUtils(qeUtils):
    def _create_input(self, sims, **kwargs):
        # This one's a bit more complicated. It's going through sims and writing the total energy and volume of each one.
        with open(self._file_path_root + self._session._input_file, "w+") as f:
            for s in sims:
                total_energy = simple_search.find_cuds_objects_by_oclass(oclass = QE.TotalEnergy, root = s, rel = QE.HAS_PART)[0].value
                volume = simple_search.find_cuds_objects_by_oclass(oclass = QE.Volume, root = s, rel = QE.HAS_PART)[0].value
                f.write(f"{volume} {total_energy}\n")
    
    def _update_cuds(self, sim):
        # Updates equilibrium volume and bulk modulus.
        with open(self._file_path_root + self._session._output_file, 'r') as file:
            lines = file.readlines()
            v0 = lines[1].split()[3]
            b0 = lines[1].split()[6][1:]
            for s in sim:
                volume_entity = s.get(oclass = QE.Cell)[0].get(oclass = QE.EquilibriumVolume)
                modulus_entity = s.get(oclass = QE.BulkModulus)
                if volume_entity:
                    volume_entity[0].value = float(v0)
                    volume_entity[0].unit = "au^3"
                else:
                    s.get(oclass = QE.Cell)[0].add(QE.EquilibriumVolume(value = v0, unit = "au^3"))
                if modulus_entity:
                    modulus_entity[0].value = float(b0)
                    volume_entity[0].unit = "kbar"
                else:
                    s.add(QE.BulkModulus(value = b0, unit = "kbar"))
        super()._update_cuds(sim)

class phUtils(qeUtils):
    def _create_input(self, sim, **kwargs):
        self.params = {
            "INPUTPH": {
                "outdir": "'.'",
                "prefix": f"'{self._session._prefix}'",
                "fildyn": f"'{self._session._prefix}.ph.dyn'"
            }
        }
        self.qpoints = []
        for point in sim.get(oclass = QE.QPoint):
            if point.calculate == True:
                self.qpoints.append(point)
        try:
            if self.params["ldisp"] != ".true." and self.params["qplot"] != ".true.":
                self.sysinfo = {
                    "": [["0 0 0"]] 
                                    # TODO: manual q point
                                    # TODO: add support for multiple q points
            }
            else:
                self.sysinfo = {}
        except:
            self.sysinfo = {}
        super()._create_input(sim, **kwargs)

    def _update_cuds(self, sim):
        sim.add(QE.PhOut(path = self._file_path_root + self._session._output_file))
        with open(self._file_path_root + self._session._output_file, 'r') as file:
            lines = file.readlines()
            beginend = []
            for i, line in enumerate(lines):
                if line.startswith(" ****"):
                    beginend.append(i)
            q_point = self.qpoints[0]
            for i in range(beginend[0]+1, beginend[1]):
                freq = float(lines[i].split()[4])
                modenum = int(lines[i].split()[2][:-1])
                unit = lines[i].split()[5][1:-1]
                for mode in q_point.get(oclass = QE.Mode):
                    if mode.number == modenum:
                        if mode.get(oclass = QE.Frequency):
                            mode.get(oclass = QE.Frequency)[0].value = freq
                            mode.get(oclass = QE.Frequency)[0].unit = unit
                        else:
                            mode.add(QE.Frequency(value = freq, unit = unit))

        super()._update_cuds(sim)

class alpha2fUtils(qeUtils):
    def _create_input(self, sim, **kwargs):
        self.params = {
            "INPUTPH": {
                "outdir": "'.'",
                "prefix": f"'{self._session._prefix}'",
                "fildyn": f"'{self._session._prefix}.ph.dyn'"
            },
            "INPUTa2F": {
                "nfreq": 500
            }
        }
        super()._create_input(sim, **kwargs)

    def _update_cuds(self, sim):

        super()._update_cuds(sim)

# class epaUtils(qeUtils):
#     def _create_input(self, sim, **kwargs):
#         with open(self._file_path_root + self._session._input_file, "w+") as file:
            