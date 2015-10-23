class IntLabel:
    def __init__(self, i):
        self._i = i

    def as_string(self):
        return str(self._i)

    def __hash__(self):
        return hash(self._i)

    def __eq__(self, other):
        return self._i == other._i

    def __repr__(self):
        return self.as_string()
        
