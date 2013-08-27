'''
Created on Aug 14, 2013

@author: caofan
'''

from exo.models import Mapping, Raw_sample, Method, Parameter_set, Genome, Parameter, Sample
from django.db import transaction, IntegrityError
from fileHandler import FileHandler
from add_params_to_db import add_params_to_db

def add_mapping_to_db():
    f = FileHandler('../../data/mapping.txt')
    for r in f:
        tokens = r.split('\t')
        rsname, gn, mn, version = tokens[0:4]
        mapped_count, unique_count, pos_unique_count, neg_unique_count, s_name = tokens[4:9]
        params = tokens[9:]
        rs = Raw_sample.objects.get(name=rsname)
        genome = Genome.objects.get(assembly_name=gn)
        method = Method.objects.get(method_name=mn, method_version=version)
        
        sid = transaction.savepoint()
        try:
            '''
            param_set = Parameter_set(method_obj = method)
            param_set.save()
            for p in params:
                t = p.split(':')
                if len(t) == 1:
                    t.append('')
                param = Parameter(param_set = param_set, param_name = t[0], param_value = t[1])
                param.save()
            '''
            param_set = add_params_to_db(params, method)
            mapping = Mapping(method_obj=method, genome_obj=genome, params=param_set, sample_obj=rs)
            mapping.save()
            sample = Sample(name = s_name, mapped_tags = mapped_count, unique_tags = unique_count, positive_unique_count = pos_unique_count, negative_unique_count = neg_unique_count, mapping_obj = mapping)
            sample.save()
            transaction.savepoint_commit(sid)
        except IntegrityError:
            print 'error'
            transaction.savepoint_rollback(sid)
    f.close()
            
        

if __name__ == '__main__':
    add_mapping_to_db()