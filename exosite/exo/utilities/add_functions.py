'''
Created on Aug 15, 2013

@author: caofan
'''

from exo.models import *
from django.db import transaction, IntegrityError
from fileHandler import FileHandler


def add_genome_to_db():
    f = FileHandler("../../data/genomes.txt")
    #sid = transaction.savepoint()
    for r in f:
        tokens = r.split('\t')
        g = Genome(assembly_name = tokens[0],
                   alternative_assembly_name = tokens[1],
                   species = tokens[2],
                   size = tokens[3])
        sid = transaction.savepoint()    #create save points before transaction for rollback if error occurs.
        try:
            g.save()
            transaction.savepoint_commit(sid)
        except IntegrityError:
            print "error"
            transaction.savepoint_rollback(sid)
    f.close()