from django.db import models


# Create your models here.

class Genome(models.Model):
    '''
    Table to store the meta data about reference genome.
    '''
    assembly_name = models.CharField(max_length=20,unique=True)
    alternative_assembly_name = models.CharField(max_length=50)
    species = models.CharField(max_length=50)
    size = models.BigIntegerField(null=False)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    def __unicode__(self):
        return self.assembly_name

class Raw_sample(models.Model):
    name = models.CharField(max_length=100,unique=True,null=False)
    tag_length = models.PositiveIntegerField(null=False)
    total_tags = models.PositiveIntegerField(null=False)
    coverage = models.FloatField(null=False)
    factor = models.CharField(max_length=20)
    antibody = models.CharField(max_length=50)
    cell_type = models.CharField(max_length=20)
    replicate_num = models.PositiveSmallIntegerField()
    
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    def __unicode__(self):
        return self.name


'''
class Mapping_method(models.Model):
    mapping_method = models.CharField(max_length=20)
    
    mapping_method_version = models.CharField(max_length=20)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    def __unicode__(self):
        return '%s_%s'%(self.mapping_method, self.mapping_method_version,)
'''

class Method(models.Model):
    method_name = models.CharField(max_length=20)
    method_type = models.CharField(max_length=20)
    method_version = models.CharField(max_length=20)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    def __unicode__(self):
        return '%s_%s'%(self.method_name, self.method_version)
    
    class Meta:
        unique_together = (('method_name','method_version'),)

class Parameter_set(models.Model):
    method_obj = models.ForeignKey(Method)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    def __unicode__(self):
        return "%s_%s_%s"%(self.method_obj.method_name, self.method_obj.method_version, self.id)
    
class Parameter(models.Model):
    param_set = models.ForeignKey(Parameter_set, null=False)
    param_name = models.CharField(max_length=50)
    param_value = models.CharField(max_length=200) #some string values can be very long)

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    class Meta:
        unique_together = (("param_set","param_name"),)
        
    def __unicode__(self):
        return '%s %s'%(self.param_set.id, self.param_name)

class Mapping(models.Model):
    method_obj = models.ForeignKey(Method)
    sample_obj = models.ForeignKey(Raw_sample)
    genome_obj = models.ForeignKey(Genome)
    params = models.ForeignKey(Parameter_set)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    def __unicode__(self):
        return "%s %s %s"%(self.sample_obj.name, self.genome_obj.assembly_name, self.method_obj)
    


'''
class Mapping_parameters(models.Model):
    mapping_obj = models.ForeignKey(Mapping)
    name = models.CharField(max_length=50)
    value = models.CharField(max_length=200)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    def __unicode__(self):
        return '%s %s'%(self.id, self.name, self.value)
'''
    
class Sample(models.Model):
    name = models.CharField(max_length=100, unique=True)
    mapped_tags = models.PositiveIntegerField(null=False)
    unique_tags = models.PositiveIntegerField(null=False)
    positive_unique_count = models.PositiveIntegerField(null=False)
    negative_unique_count = models.PositiveIntegerField(null=False)
    mapping_obj = models.OneToOneField(Mapping)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    def __unicode__(self):
        return '%s %s'%(self.name, self.mapping_obj.method_obj.method_name)


class Call_peak(models.Model):
    sample_obj = models.ForeignKey(Sample)
    method_obj = models.ForeignKey(Method)
    params = models.ForeignKey(Parameter_set)
    
    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

'''
class Peak_call_parameter(models.Model):
    run = models.ForeignKey(Call_peak
    name = models.CharField(max_length=50)
    value = models.CharField(max_length=200) #some string values can be very long)

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    class Meta:
        unique_together = (("run","name"),)
'''



class Peak(models.Model):
    run = models.ForeignKey(Call_peak)
    chrom = models.CharField(max_length=40)
    start = models.PositiveIntegerField( null=False ) #0-based
    end = models.PositiveIntegerField( null=False ) #1-based
    size = models.PositiveIntegerField(null=False)
    strand = models.CharField(max_length=1, choices=(('+','Positive'),('-', 'Negative'),('.','Unspecified'),))
    tag_count = models.PositiveIntegerField(null=False)
    summit_pos_pileup = models.IntegerField(null=False) #Origin summit location by peak calling method
    summit_val_pileup = models.FloatField(null=False)
    summit_pos_5 = models.IntegerField()
    summit_val_5 = models.FloatField()
    p_value = models.FloatField()
    q_value = models.FloatField()
    fold = models.FloatField()

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    def __unicode__(self):
        return '%s_%s_%s_%s'%(self.id, self.chrom, self.start, self.end)



class Pairing(models.Model):
    run = models.ForeignKey(Call_peak)
    method_obj = models.ForeignKey(Method)
    params = models.ForeignKey(Parameter_set)

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    def __unicode__(self):
        return self.id

'''
class Pairing_parameter(models.Model):
    pairing_obj = models.ForeignKey(Pairing)
    name = models.CharField(max_length=50)
    value = models.CharField(max_length=200)

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    class Meta:
        unique_together = (("pairing_obj","name"),)
'''




class Peak_pair(models.Model):
    '''
    Parameters for pairing peaks.
    '''
    pairing = models.ForeignKey(Pairing)
    peak1 = models.ForeignKey(Peak, related_name='+')
    peak2 = models.ForeignKey(Peak, related_name='+')
    shift_distance = models.PositiveIntegerField(null=False)
    total_score = models.FloatField(null=False)
    average_height = models.FloatField()
    summit_height = models.FloatField()

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)

    class Meta:
        unique_together = (("pairing","peak1","peak2"),)


class Chromosome(models.Model):
    genome_obj = models.ForeignKey(Genome)
    name = models.CharField(max_length=30)
    length = models.PositiveIntegerField(null=False)

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    class Meta:
        unique_together = (("genome_obj","name"),)
        
    def __unicode__(self):
        return '%s_%s'%(self.genome_obj.assembly_name, self.name,)


class Gene(models.Model):
    '''
    The table to store all genes.
    '''
    chromosome_obj = models.ForeignKey(Chromosome)
    identifier = models.CharField(max_length=50)
    name = models.CharField(max_length=30)
    start_pos = models.PositiveIntegerField(null=False) #0-based
    end_pos = models.PositiveIntegerField(null=False) #1-based
    source = models.CharField(max_length=100) #where is the gene annotation obtained. Should include the version of the source.

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    class Meta:
        unique_together = (("chromosome_obj","identifier"),("chromosome_obj","name"),)

class Peak_pair_annotation(models.Model):
    peak_pair = models.ForeignKey(Peak_pair)
    gene = models.ForeignKey(Gene)
    distance = models.PositiveIntegerField(null=False)
    type = models.CharField(max_length=30)

    created_at = models.DateTimeField(auto_now_add = True)
    updated_at = models.DateTimeField(auto_now = True)
    
    class Meta:
        unique_together = (("peak_pair","gene"),)
