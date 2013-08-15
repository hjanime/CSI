'''
Created on Aug 14, 2013

@author: caofan
'''

from exo.models import Mapping, Raw_sample, Method, Parameter_set, Genome, Parameter
from django.db import transaction, IntegrityError
from fileHandler import FileHandler

def add_mapping_to_db():
    f = FileHandler('../../data/mapping.txt')
    for r in f:
        tokens = r.split('\t')
        rsname, gn, mn, version = tokens[0:4]
        params = tokens[4:]
        rs = Raw_sample.objects.get(name=rsname)
        genome = Genome.objects.get(assembly_name=gn)
        method = Method.objects.get(method_name=mn, method_version=version)
        
        sid = transaction.savepoint()
        try:
            param_set = Parameter_set(method_obj = method)
            param_set.save()
            for p in params:
                t = p.split(':')
                if len(t) == 1:
                    t.append('')
                param = Parameter(param_set = param_set, param_name = t[0], param_value = t[1])
                param.save()
            mapping = Mapping(method_obj=method, genome_obj=genome, params=param_set, sample_obj=rs)
            mapping.save()
            transaction.savepoint_commit(sid)
        except IntegrityError:
            print 'error'
            transaction.savepoint_rollback(sid)
    f.close()
            
        

if __name__ == '__main__':
    add_mapping_to_db()