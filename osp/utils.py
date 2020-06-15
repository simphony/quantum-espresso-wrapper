
class qeInputScript:

    SystemSections = ('&CONTROL', '&SYSTEM', '&ELECTRONS')

    ParameterSections = ('K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS')

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
                output[section].append(line)
        return output