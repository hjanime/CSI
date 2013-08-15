# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Genome'
        db.create_table(u'exo_genome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('assembly_name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=20)),
            ('alternative_assembly_name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('species', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('size', self.gf('django.db.models.fields.BigIntegerField')()),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Genome'])

        # Adding model 'Raw_sample'
        db.create_table(u'exo_raw_sample', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=100)),
            ('tag_length', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('total_tags', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('coverage', self.gf('django.db.models.fields.FloatField')()),
            ('factor', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('antibody', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('cell_type', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('replicate_num', self.gf('django.db.models.fields.PositiveSmallIntegerField')()),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Raw_sample'])

        # Adding model 'Method'
        db.create_table(u'exo_method', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('method_name', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('method_type', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('method_version', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Method'])

        # Adding unique constraint on 'Method', fields ['method_name', 'method_version']
        db.create_unique(u'exo_method', ['method_name', 'method_version'])

        # Adding model 'Parameter_set'
        db.create_table(u'exo_parameter_set', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('method_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Method'])),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Parameter_set'])

        # Adding model 'Parameter'
        db.create_table(u'exo_parameter', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('param_set', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Parameter_set'])),
            ('param_name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('param_value', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Parameter'])

        # Adding unique constraint on 'Parameter', fields ['param_set', 'param_name']
        db.create_unique(u'exo_parameter', ['param_set_id', 'param_name'])

        # Adding model 'Mapping'
        db.create_table(u'exo_mapping', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('method_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Method'])),
            ('sample_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Raw_sample'])),
            ('genome_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Genome'])),
            ('params', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Parameter_set'])),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Mapping'])

        # Adding model 'Sample'
        db.create_table(u'exo_sample', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('mapped_tags', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('unique_tags', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('positive_unique_count', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('negative_unique_count', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('mapping_obj', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['exo.Mapping'], unique=True)),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Sample'])

        # Adding model 'Call_peak'
        db.create_table(u'exo_call_peak', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Sample'])),
            ('method_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Method'])),
            ('params', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Parameter_set'])),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Call_peak'])

        # Adding model 'Peak'
        db.create_table(u'exo_peak', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('run', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Call_peak'])),
            ('chrom', self.gf('django.db.models.fields.CharField')(max_length=40)),
            ('start', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('end', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('size', self.gf('django.db.models.fields.PositiveSmallIntegerField')()),
            ('strand', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('tag_count', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('summit_pos_pileup', self.gf('django.db.models.fields.PositiveSmallIntegerField')()),
            ('summit_val_pileup', self.gf('django.db.models.fields.FloatField')()),
            ('summit_pos_5', self.gf('django.db.models.fields.PositiveSmallIntegerField')()),
            ('summit_val_5', self.gf('django.db.models.fields.FloatField')()),
            ('p_value', self.gf('django.db.models.fields.FloatField')()),
            ('q_value', self.gf('django.db.models.fields.FloatField')()),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Peak'])

        # Adding model 'Pairing'
        db.create_table(u'exo_pairing', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('run', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Call_peak'])),
            ('method_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Method'])),
            ('params', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Parameter_set'])),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Pairing'])

        # Adding model 'Peak_pair'
        db.create_table(u'exo_peak_pair', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('pairing', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Pairing'])),
            ('peak1', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['exo.Peak'])),
            ('peak2', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['exo.Peak'])),
            ('shift_distance', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('total_score', self.gf('django.db.models.fields.FloatField')()),
            ('average_height', self.gf('django.db.models.fields.FloatField')()),
            ('summit_height', self.gf('django.db.models.fields.FloatField')()),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Peak_pair'])

        # Adding unique constraint on 'Peak_pair', fields ['pairing', 'peak1', 'peak2']
        db.create_unique(u'exo_peak_pair', ['pairing_id', 'peak1_id', 'peak2_id'])

        # Adding model 'Chromosome'
        db.create_table(u'exo_chromosome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('genome_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Genome'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Chromosome'])

        # Adding unique constraint on 'Chromosome', fields ['genome_obj', 'name']
        db.create_unique(u'exo_chromosome', ['genome_obj_id', 'name'])

        # Adding model 'Gene'
        db.create_table(u'exo_gene', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('chromosome_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Chromosome'])),
            ('identifier', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('start_pos', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('end_pos', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('source', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Gene'])

        # Adding unique constraint on 'Gene', fields ['chromosome_obj', 'identifier']
        db.create_unique(u'exo_gene', ['chromosome_obj_id', 'identifier'])

        # Adding unique constraint on 'Gene', fields ['chromosome_obj', 'name']
        db.create_unique(u'exo_gene', ['chromosome_obj_id', 'name'])

        # Adding model 'Peak_pair_annotation'
        db.create_table(u'exo_peak_pair_annotation', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('peak_pair', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Peak_pair'])),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Gene'])),
            ('distance', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('created_at', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('updated_at', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'exo', ['Peak_pair_annotation'])

        # Adding unique constraint on 'Peak_pair_annotation', fields ['peak_pair', 'gene']
        db.create_unique(u'exo_peak_pair_annotation', ['peak_pair_id', 'gene_id'])


    def backwards(self, orm):
        # Removing unique constraint on 'Peak_pair_annotation', fields ['peak_pair', 'gene']
        db.delete_unique(u'exo_peak_pair_annotation', ['peak_pair_id', 'gene_id'])

        # Removing unique constraint on 'Gene', fields ['chromosome_obj', 'name']
        db.delete_unique(u'exo_gene', ['chromosome_obj_id', 'name'])

        # Removing unique constraint on 'Gene', fields ['chromosome_obj', 'identifier']
        db.delete_unique(u'exo_gene', ['chromosome_obj_id', 'identifier'])

        # Removing unique constraint on 'Chromosome', fields ['genome_obj', 'name']
        db.delete_unique(u'exo_chromosome', ['genome_obj_id', 'name'])

        # Removing unique constraint on 'Peak_pair', fields ['pairing', 'peak1', 'peak2']
        db.delete_unique(u'exo_peak_pair', ['pairing_id', 'peak1_id', 'peak2_id'])

        # Removing unique constraint on 'Parameter', fields ['param_set', 'param_name']
        db.delete_unique(u'exo_parameter', ['param_set_id', 'param_name'])

        # Removing unique constraint on 'Method', fields ['method_name', 'method_version']
        db.delete_unique(u'exo_method', ['method_name', 'method_version'])

        # Deleting model 'Genome'
        db.delete_table(u'exo_genome')

        # Deleting model 'Raw_sample'
        db.delete_table(u'exo_raw_sample')

        # Deleting model 'Method'
        db.delete_table(u'exo_method')

        # Deleting model 'Parameter_set'
        db.delete_table(u'exo_parameter_set')

        # Deleting model 'Parameter'
        db.delete_table(u'exo_parameter')

        # Deleting model 'Mapping'
        db.delete_table(u'exo_mapping')

        # Deleting model 'Sample'
        db.delete_table(u'exo_sample')

        # Deleting model 'Call_peak'
        db.delete_table(u'exo_call_peak')

        # Deleting model 'Peak'
        db.delete_table(u'exo_peak')

        # Deleting model 'Pairing'
        db.delete_table(u'exo_pairing')

        # Deleting model 'Peak_pair'
        db.delete_table(u'exo_peak_pair')

        # Deleting model 'Chromosome'
        db.delete_table(u'exo_chromosome')

        # Deleting model 'Gene'
        db.delete_table(u'exo_gene')

        # Deleting model 'Peak_pair_annotation'
        db.delete_table(u'exo_peak_pair_annotation')


    models = {
        u'exo.call_peak': {
            'Meta': {'object_name': 'Call_peak'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Method']"}),
            'params': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Parameter_set']"}),
            'sample_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Sample']"}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.chromosome': {
            'Meta': {'unique_together': "(('genome_obj', 'name'),)", 'object_name': 'Chromosome'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'genome_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Genome']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.gene': {
            'Meta': {'unique_together': "(('chromosome_obj', 'identifier'), ('chromosome_obj', 'name'))", 'object_name': 'Gene'},
            'chromosome_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Chromosome']"}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'end_pos': ('django.db.models.fields.PositiveIntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'identifier': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'start_pos': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.genome': {
            'Meta': {'object_name': 'Genome'},
            'alternative_assembly_name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'assembly_name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '20'}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'size': ('django.db.models.fields.BigIntegerField', [], {}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.mapping': {
            'Meta': {'object_name': 'Mapping'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'genome_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Genome']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Method']"}),
            'params': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Parameter_set']"}),
            'sample_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Raw_sample']"}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.method': {
            'Meta': {'unique_together': "(('method_name', 'method_version'),)", 'object_name': 'Method'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_name': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'method_type': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'method_version': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.pairing': {
            'Meta': {'object_name': 'Pairing'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Method']"}),
            'params': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Parameter_set']"}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Call_peak']"}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.parameter': {
            'Meta': {'unique_together': "(('param_set', 'param_name'),)", 'object_name': 'Parameter'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'param_name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'param_set': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Parameter_set']"}),
            'param_value': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.parameter_set': {
            'Meta': {'object_name': 'Parameter_set'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Method']"}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.peak': {
            'Meta': {'object_name': 'Peak'},
            'chrom': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'end': ('django.db.models.fields.PositiveIntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'p_value': ('django.db.models.fields.FloatField', [], {}),
            'q_value': ('django.db.models.fields.FloatField', [], {}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Call_peak']"}),
            'size': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'start': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'summit_pos_5': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'summit_pos_pileup': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'summit_val_5': ('django.db.models.fields.FloatField', [], {}),
            'summit_val_pileup': ('django.db.models.fields.FloatField', [], {}),
            'tag_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.peak_pair': {
            'Meta': {'unique_together': "(('pairing', 'peak1', 'peak2'),)", 'object_name': 'Peak_pair'},
            'average_height': ('django.db.models.fields.FloatField', [], {}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'pairing': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Pairing']"}),
            'peak1': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['exo.Peak']"}),
            'peak2': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['exo.Peak']"}),
            'shift_distance': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'summit_height': ('django.db.models.fields.FloatField', [], {}),
            'total_score': ('django.db.models.fields.FloatField', [], {}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.peak_pair_annotation': {
            'Meta': {'unique_together': "(('peak_pair', 'gene'),)", 'object_name': 'Peak_pair_annotation'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'distance': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'peak_pair': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Peak_pair']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.raw_sample': {
            'Meta': {'object_name': 'Raw_sample'},
            'antibody': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'cell_type': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'coverage': ('django.db.models.fields.FloatField', [], {}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'factor': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '100'}),
            'replicate_num': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'tag_length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'total_tags': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'exo.sample': {
            'Meta': {'object_name': 'Sample'},
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'mapped_tags': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'mapping_obj': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['exo.Mapping']", 'unique': 'True'}),
            'negative_unique_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'positive_unique_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'unique_tags': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['exo']