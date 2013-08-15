'''
Created on Aug 14, 2013

@author: caofan
'''

class FileHandler:
    def __init__(self, filename):
        self.filename = filename
        self.handler = open(filename)
    
    
    def __iter__(self):
        return self
    
    def next(self):
        r = self.handler.next()
        r = r.strip()
        if r != '' and r[0] != '#':
            return r
        else:
            return self.next()
    
    
    def close(self):
        self.handler.close()
    