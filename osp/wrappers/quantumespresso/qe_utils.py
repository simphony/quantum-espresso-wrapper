from osp.core import QE

class qeUtils():
    """Utilities for reading and writing .in and .out files
    """
    SystemSections = ('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS')


    def __init__(self, session):
        self._session = session
    # def __init__(self):
    #     return None

    def _read_input(self, file_path):
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
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("!"):
                    energy = float(line.split()[4])
        try: 
            return energy
        except:
            print("Error: output file did not contain total energy!")

    def _create_input(self, file_path, simulation):
        CONTROL = ["&CONTROL", 
        "   calculation = 'scf',", 
        "   outdir = '.',", 
        "   prefix = 'basic',", 
        "   pseudo_dir = '.',", 
        "   verbosity = 'low',", 
        "   tprnfor = .true.,", 
        "   tstress=.true.", 
        "/"]
        SYSTEM = ["&SYSTEM",
        "   ibrav = 0",
        "   A = 5.43070",
        "   nat = 2",
        "   ntyp = 1",
        "   ecutwfc=50,",
        "   ecutrho=200,",
        "   input_dft='pbe',",
        "   occupations='smearing',",
        "   smearing='mv',",
        "   degauss=0.005d0,",
        "/"]
        ELECTRONS = ["&ELECTRONS",
        "   conv_thr=1d-08,",
	    "   mixing_beta=0.7d0,",
        "/"]
        with open(file_path, 'w+') as file:
            for line in CONTROL:
                file.write(line + '\n')
            for line in SYSTEM:
                file.write(line + '\n')
            for line in ELECTRONS:
                file.write(line + '\n')
            file.write("CELL_PARAMETERS {alat}\n")
            for cell in simulation.get(oclass = QE.CELL):
                for cell_params in cell.get(oclass = QE.CELL_PARAMS):
                    for param in cell_params.iter():
                        file.write("  " + str(param.vector[0]) + " " + str(param.vector[1]) + " " + str(param.vector[2]) + "\n")

            file.write("ATOMIC_SPECIES \n")
            for element in simulation.get(oclass = QE.ELEMENT):
                file.write(" " + element.name)
                for mass in element.get(oclass = QE.MASS):
                    file.write(" " + str(mass.value))
                for pseudo in element.get(oclass = QE.PSEUDOPOTENTIAL):
                    file.write(" " + pseudo.name + "\n")
            file.write("ATOMIC_POSITIONS {crystal} \n")
            for element in simulation.get(oclass = QE.ELEMENT):
                for atom in element.get(oclass = QE.ATOM):
                    for position in atom.get(oclass = QE.POSITION):
                        file.write(" " + element.name + " " + str(position.vector[0]) + " " + str(position.vector[1]) + " " + str(position.vector[2]) + "\n")
            file.write("K_POINTS {automatic}\n  7 7 7 0 0 0")


            

        
