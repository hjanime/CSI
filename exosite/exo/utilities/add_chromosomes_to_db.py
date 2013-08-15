'''
Created on Aug 13, 2013

@author: caofan
'''

from exo.models import Genome, Chromosome

from django.db import transaction, IntegrityError
from fileHandler import FileHandler

def add_chromosomes_to_db():
    genomes = Genome.objects.all()
    if len(genomes) <= 0:
        print "No genomes. Abort."
    else:
        for g in genomes:
            f = FileHandler('../../data/%s.chrom.sizes'%g.assembly_name)
            #sid = transaction.savepoint()
            for r in f:
                cname, csize = r.split()
                chrom = Chromosome(name=cname, length=csize, genome_obj=g)
                sid = transaction.savepoint()
                try:
                    chrom.save()
                    transaction.savepoint_commit(sid)
                except IntegrityError:
                    transaction.savepoint_rollback(sid)
                    print 'error'
            
            f.close()


if __name__ == '__main__':
    add_chromosomes_to_db()
                