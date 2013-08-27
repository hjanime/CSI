'''
Created on 20 Aug, 2013

@author: caofan
'''
from exo.models import Parameter, Parameter_set

def add_params_to_db(params, method):
    param_set = Parameter_set(method_obj = method)
    param_set.save()
    for p in params:
        t = p.split(':')
        if len(t) == 1:
            t.append('')
        param = Parameter(param_set = param_set, param_name = t[0], param_value = t[1])
        param.save()
    return param_set
