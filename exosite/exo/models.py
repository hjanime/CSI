from django.db import models

# Create your models here.
class Sample(models.Model):
    name = models.CharField(max_length=100,unique=True,null=False)
    tag_length = models.PositiveIntegerField(null=False)
    total_tags = models.PositiveIntegerField(null=False)
    mapped_tags = models.PositiveIntegerField(null=False)
    unique_tags = models.PositiveIntegerField(null=False)
    factor = models.CharField(max_length=20)
    antibody = models.CharField(max_length=50)
    replicate_num = models.PositiveSmallInteger()
    cell = models.CharField(max_length=20)
    reference_genome = models.CharField(max_length=50)


class Sample_run(models.Model):
    sample = models.ForeignKey(Sample)
    shift_size = models.PositiveSmallInteger(null=False)
    tag_size = models.PositiveSmallInteger(null=False)
    method = models.CharField(max_length=15)
    version = models.CharField(max_length=15)
    p_value_thresh = models.FloatField()
    q_value_thresh = models.FloatField()
    rpm_thresh = models.FloatField()
    rpkm_thresh = models.FloatField()
    pileup_thresh = models.FloatField()
    fdr_thresh = models.FloatField()


class Peak(models.Model):
    run = models.ForeignKey(Sample_run)
    chrom = models.CharField(max_length=40)
    start = models.PositiveInteger( null=False ) #0-based
    end = models.PositiveInteger( null=False ) #1-based
    size = models.PositiveSmallInteger(null=False)
    strand = models.CharField(max_length=10) #add choices
    tag_count = models.PositiveInteger(null=False)
    summit_pos_pileup = models.PositiveSmallInteger(null=False) #Origin summit location by peak calling method
    summit_val_pileup = models.FloatField(null=False)
    summit_pos_5 = models.PositiveSmallInteger()
    summit_val_5 = models.FloatField()
    p_value = models.FloatField()
    q_value = models.FloatField()


class Pair(models.Model):
    '''
    Parameters for pairing peaks.
    '''
    

