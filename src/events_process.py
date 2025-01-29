import pandas as pd
from events import Events

class EventsProcess(Events):

    def __init__(self, event):
        super().__init__(event)

        ## need to manually set the fields attribute
        self.stimtype = None 
        self.stimparam = None 
        self.stimDescription = None 
        self.OEStartTime = None
        self.OEEndTime = None 
        self.freqStart = None
        self.freqEnd = None
        self.firstFrame = None 
        self.lastFrame = None
        self.frameRate = None
        self.spikeTimes = []
        self.pupilDiameter = []
        self.frameNums = []
        self.frameTimes = []
        self.timeWindow = None 
        self.frameWindow = None 

    def toDataFrame(self):
        data = {field: [] for field in self.fields}
        for field in self.fields:
            data[field].append(getattr(self, field))
        return pd.DataFrame(data)
    

    def __repr__(self):
        return f"EventsProcess({', '.join([f'{field}={getattr(self, field)}' for field in self.fields])}"
