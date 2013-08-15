'''
Created on Aug 14, 2013

@author: caofan
'''

from exo.models import Raw_sample
from django.db import transaction, DatabaseError
from fileHandler import FileHandler


def add_raw_to_db():
    f = FileHandler("../../data/raw_samples.txt")
    #sid = transaction.savepoint()
    for r in f:
        tokens = r.split('\t')
        rs = Raw_sample(name = tokens[0],
                        tag_length = tokens[1],
                        total_tags = tokens[2],
                        coverage = tokens[3],
                        factor = tokens[4],
                        antibody = tokens[5],
                        cell_type = tokens[6],
                        replicate_num = tokens[7])
        sid = transaction.savepoint()
        try:
            rs.save()
            transaction.savepoint_commit(sid)
        except DatabaseError, msg:
            print msg
            transaction.savepoint_rollback(sid)
    f.close()
                
if __name__ == '__main__':
    add_raw_to_db()