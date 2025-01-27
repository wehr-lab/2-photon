from utils.extract_vals import extract_value


class Cell:
    def __init__(self, sortedUnit):
        if not hasattr(sortedUnit, 'dtype') or not sortedUnit.dtype.names:
            raise ValueError('Input should be a numpy record array with named fields.')
        
        self.fields = sortedUnit.dtype.names

        ## set attributes for each field in the event

        for field in self.fields:
            value = extract_value(sortedUnit[field], field)
            setattr(self, field, value)

        def __repr__(self):
            return f'Cell {self.cellnum}'