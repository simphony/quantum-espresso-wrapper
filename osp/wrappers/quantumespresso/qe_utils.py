
class qeInputUtils():
    """Utilities for reading and writing .in files
    """
    SystemSections = ('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS')


    def __init__(self, filepath):


    def _read_file(self, file_path):

        with open(file_path, "r") as f:


    def _split_into_sections(self, lines):
        output = {}
        for line in lines:
            if line.startswith(self.SystemSections):
                section = line
                output[section] = []
            else:
                if line != "/":
                    output[section].append(line)
                
    def _write_file(self, content, file_path):
        


class qeInputScript:

    ParameterSections = ('&CONTROL', '&SYSTEM', '&ELECTRONS')

    SystemSections = ('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS')

    def __init__(self, filename):
        self._filename = filename
        self._content = None

    def parse(self):
        text = self._read_clean()
        self._content = self._split_into_sections(text)

    def _read_clean(self):
        content = []
        with open(self._filename, 'r') as f:
            for line in f.readlines():
                line = line.partition('#')[0].strip()
                if line:
                    content.append(line)
        return content
    
    def _split_into_sections(self, lines):
        output = {}
        section = 'Header'
        output[section] = []
        for line in lines:
            if line.startswith(self.ParameterSections) or line.startswith(self.SystemSections):
                section = line
                output[section] = []
            else:
                if line != "/":
                    output[section].append(line)
        return output


