'''
Created on Aug 14, 2013

@author: caofan
'''

from exo.models import Method
from fileHandler import FileHandler
from django.db import transaction, DatabaseError

def add_method_to_db():
    f = FileHandler("../../data/methods.txt")
    sid = transaction.savepoint()
    for r in f:
        tokens = r.split('\t')
        m = Method(method_name=tokens[0],
                   method_type=tokens[1],
                   method_version=tokens[2])
        sid = transaction.savepoint()
        try:
            m.save()
            transaction.savepoint_commit(sid)
        except DatabaseError:
            transaction.savepoint_rollback(sid)
    f.close()
        

if __name__=='__main__':
    add_method_to_db()