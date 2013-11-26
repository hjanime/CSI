#code from http://www.moosechips.com/2011/01/python-chip-seq-bed-file-reader/
#/usr/bin/python

import csv
class CommentedFileReader:
    """
    Helper class for file reading.
    Skips lines starting with '#'

    tsv_file = csv.reader(CommentedFileReader("inputfile.txt"),
                      delimiter='\t')
    for row in tsv_file:
        print row[2] # prints column 3 of each line
    """
    def __init__(self, f, commentstring="#"):
        self.f = open(f, 'rU')
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring) or line.startswith('track'):
            line = self.f.next()
        return line
    def __iter__(self):
        return self

csv.register_dialect('bed', delimiter = '\t',quoting = csv.QUOTE_NONE,skipinitialspace = True)

class BEDReader:
    """
    Read BED files into a DictReader.
    Usage: BEDReader( filename, filetype ).
    filename: The name of the file to read.
    filetype: As there are currently multiple extensions
              to the original Bed format, we included them.
              The field can be: 'bed12', 'bed3','bed6', 'bedGraph',
              'narrowPeak', 'broadPeak'.

    Use BEDReader.fieldnames to see the field names.

    Example:
        bed = BEDReader("file.bed", "bed12")
        for line in bed:
            # print the chromStart
            print(line['chromStart'])
    """
    FIELDS = ['chrom', 'chromStart', 'chromEnd',
              'name', 'score', 'strand',
              'thickStart', 'thickEnd',
              'itemRgb',
              'blockCount', 'blockSizes', 'blockStarts', 'value', 'signalValue', 'pValue', 'qValue', 'peak']

    def __init__(self, filename, ftype):
        ftype = ftype.strip().lower()
        if ftype == 'bed12':
            fieldNames = self.FIELDS[:12]
        elif ftype == 'bed3':
            fieldNames = self.FIELDS[:3]
        elif ftype == 'bed6':
            fieldNames = self.FIELDS[:6]
        elif ftype == 'bedgraph':
            fieldNames = self.FIELDS[:3]+self.FIELDS[-5:-4]
        elif ftype == 'narrowpeak':
            fieldNames=self.FIELDS[:6]+self.FIELDS[-4:]
        elif ftype == 'broadpeak':
            fieldNames=self.FIELDS[:6]+self.FIELDS[-4:-1]

        self.fieldnames = fieldNames
        self.reader = csv.reader(CommentedFileReader(filename), dialect='bed')

    def __iter__(self):
        return self

    def getFields(self):
        return self.fieldnames

    def next(self):
        row = self.reader.next()


        ret = []
        for r in row:
            if r.isdigit():
                ret.append(int(r))
            else:
                tmp = r
                try:
                    tmp = float(r)
                except ValueError, TypeError:
                    tmp = r
                ret.append(tmp)
        return ret

    def tell(self):
        return self.reader.tell()

    def seek(self, pos):
        return self.reader.seek(pos)




class BEDWriter(csv.DictWriter):
    FIELDS = ('chrom', 'chromStart', 'chromEnd',
              'name', 'score', 'strand',
              'thickStart', 'thickEnd',
              'itemRgb',
              'blockCount', 'blockSizes', 'blockStarts')

    def __init__(self, filename):
        csv.DictWriter.__init__(self, open(filename,'wb'), dialect='bed',fieldnames=self.FIELDS,extrasaction='ignore',quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    bed = BEDReader("4/rawbin.refgene")
    for line in bed:
        print(line['chromStart'])
