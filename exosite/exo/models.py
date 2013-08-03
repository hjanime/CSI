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

    def __unicode__(self):
        return self.assembly_name

class Raw_sample(models.Model):
    name = models.CharField(max_length=100,unique=True,null=False)
    tag_length = models.PositiveIntegerField(null=False)
    total_tags = models.PositiveIntegerField(null=False)
    factor = models.CharField(max_length=20)
    antibody = models.CharField(max_length=50)
    replicate_num = models.PositiveSmallIntegerField()
    cell_type = models.CharField(max_length=20)

    def __unicode__(self):
        return self.name



class Mapping_method(models.Model):
    mapping_method = models.CharField(max_length=20)
    mapping_method_version = models.CharField(max_length=20)

    def __unicode__(self):
        return '%s_%s'%(self.mapping_method, self.mapping_method_version,)

class Mapping(models.Model):
    method_obj = models.ForeignKey(Mapping_method)
    sample_obj = models.ForeignKey(Raw_sample)

class Mapping_parameters(models.Model):
    mapping_obj = models.ForeignKey(Mapping)
    name = models.CharField(max_length=50)
    value = models.CharField(max_length=200)

    def __unicode__(self):
        return '%s %s'%(self.id, self.name, self.value)

class Sample(models.Model):
    mapped_tags = models.PositiveIntegerField(null=False)
    unique_tags = models.PositiveIntegerField(null=False)
    positive_unique_count = models.PositiveIntegerField(null=False)
    negative_unique_count = models.PositiveIntegerField(null=False)
    reference_genome = models.ForeignKey(Genome)
    mapping_obj = models.ForeignKey(Mapping)


class Sample_run(models.Model):
    sample_obj = models.ForeignKey(Sample)
    method = models.CharField(max_length=20)
    version = models.CharField(max_length=20)


class Peak_call_parameter(models.Model):
    run = models.ForeignKey(Sample_run)
    name = models.CharField(max_length=50)
    value = models.CharField(max_length=200) #some string values can be very long)

    class Meta:
        unique_together = (("run","name"),)


class Pairing(models.Model):
    run = models.ForeignKey(Sample_run)
    method = models.CharField(max_length=20)
    version = models.CharField(max_length=20)



class Pairing_parameter(models.Model):
    pairing_obj = models.ForeignKey(Pairing)
    name = models.CharField(max_length=50)
    value = models.CharField(max_length=200)

    class Meta:
        unique_together = (("pairing_obj","name"),)


class Peak(models.Model):
    run = models.ForeignKey(Sample_run)
    chrom = models.CharField(max_length=40)
    start = models.PositiveIntegerField( null=False ) #0-based
    end = models.PositiveIntegerField( null=False ) #1-based
    size = models.PositiveSmallIntegerField(null=False)
    strand = models.CharField(max_length=1, choices=(('+','Positive'),('-', 'Negative'),('.','Unspecified'),))
    tag_count = models.PositiveIntegerField(null=False)
    summit_pos_pileup = models.PositiveSmallIntegerField(null=False) #Origin summit location by peak calling method
    summit_val_pileup = models.FloatField(null=False)
    summit_pos_5 = models.PositiveSmallIntegerField()
    summit_val_5 = models.FloatField()
    p_value = models.FloatField()
    q_value = models.FloatField()




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


    class Meta:
        unique_together = (("pairing","peak1","peak2"),)


class Chromosome(models.Model):
    genome_obj = models.ForeignKey(Genome)
    name = models.CharField(max_length=30)
    length = models.PositiveIntegerField(null=False)

    class Meta:
        unique_together = (("genome_obj","name"),)


class Gene(models.Model):
    '''
    The table to store all genes.
    '''
    genome_obj = models.ForeignKey(Genome)
    identifier = models.CharField(max_length=50)
    name = models.CharField(max_length=30)
    start_pos = models.PositiveIntegerField(null=False) #0-based
    end_pos = models.PositiveIntegerField(null=False) #1-based
    source = models.CharField(max_length=100) #where is the gene annotation obtained. Should include the version of the source.

    class Meta:
        unique_together = (("genome_obj","identifier"),("genome_obj","name"),)

class Peak_pair_annotation(models.Model):
    peak_pair = models.ForeignKey(Peak_pair)
    gene = models.ForeignKey(Gene)
    distance = models.PositiveIntegerField(null=False)
    type = models.CharField(max_length=30)

    class Meta:
        unique_together = (("peak_pair","gene"),)
