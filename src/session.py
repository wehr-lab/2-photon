import pandas as pd
import numpy as np
from events import Events

class Session:

    def __init__(self):
        self.events = []

    def add_event(self, event):
        if isinstance(event, Events):
            self.events.append(event)
        else:
            raise ValueError('Input should be an instance of the Events class.')
    
    def get_event(self, index):
        for event in self.events:
            if event.index == index:
                return event
            else:
                raise ValueError('Event not found.')    

    def __len__(self):
        return len(self.events)
    
    def toDataFrame(self):
        data = {field: [] for field in self.events[0].fields}
        for event in self.events:
            for field in event.fields:
                data[field].append(getattr(event, field))
        return pd.DataFrame(data)         
         

    def __repr__(self):
        return f'Session with {len(self.events)} events'
