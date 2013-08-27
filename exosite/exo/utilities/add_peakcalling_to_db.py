'''
Created on 20 Aug, 2013

@author: caofan
'''

from exo.models import  Method, Sample, Call_peak, Peak
from django.db import transaction, IntegrityError
from fileHandler import FileHandler
from add_params_to_db import add_params_to_db

def add_peakcalling_to_db():
    f = FileHandler("../../data/peakcalling.txt")
    
    
    for r in f:
        tokens = r.split('\t')
        s_name, method_name, method_version = tokens[0:3]
        peak_file = tokens[3]
        params = tokens[4:]
        print "a"
        sample = Sample.objects.get(name=s_name)
        method = Method.objects.get(method_name = method_name, method_version = method_version)
        print 'b'
        '''
        run = models.ForeignKey(Call_peak)
        chrom = models.CharField(max_length=40)
        start = models.PositiveIntegerField(null=False) #0-based
        end = models.PositiveIntegerField(null=False) #1-based
        size = models.PositiveIntegerField(null=False)
        strand = models.CharField(max_length=1, choices=(('+', 'Positive'), ('-', 'Negative'), ('.', 'Unspecified')))
        tag_count = models.FloatField(null=False)
        summit_pos_pileup = models.PositiveSmallIntegerField(null=False) #Origin summit location by peak calling method
        summit_val_pileup = models.FloatField(null=False)
        summit_pos_5 = models.PositiveSmallIntegerField()
        summit_val_5 = models.FloatField()
        p_value = models.FloatField()
        q_value = models.FloatField()
        fold = models.FloatField()
        '''
        
        
        sid = transaction.savepoint()
        
        try:
            param_set = add_params_to_db(params, method)
            callpeak = Call_peak(sample_obj = sample, method_obj = method, params = param_set)
            callpeak.save()
            sid = transaction.savepoint()
            ph = FileHandler(peak_file)
            for p in ph:
                ptokens = p.split('\t')
                chrom, start, end, length, strand, abs_summit, pileup, pileup_5, p_val, fold, q_val, _, tc, rpkm = ptokens
                start = int(start)
                abs_summit = int(abs_summit)
                tc = float(tc)
                #print p
                newp = Peak(run=callpeak, chrom=chrom, start=start-1, end=end,
                            size=length, strand=strand, tag_count=tc, summit_pos_pileup=abs_summit-1-start,
                            summit_val_pileup = pileup, summit_pos_5 = abs_summit-1-start, summit_val_5=pileup_5,
                            p_value = p_val, q_value = q_val, fold = fold)
                newp.save()
            transaction.savepoint_commit(sid)
            ph.close()
        except IntegrityError, msg:
            print msg
            transaction.savepoint_rollback(sid)
        

if __name__=='__main__':
    add_peakcalling_to_db()