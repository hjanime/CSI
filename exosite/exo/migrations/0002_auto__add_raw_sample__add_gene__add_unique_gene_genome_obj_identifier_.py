# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Raw_sample'
        db.create_table(u'exo_raw_sample', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=100)),
            ('tag_length', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('total_tags', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('factor', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('antibody', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('replicate_num', self.gf('django.db.models.fields.PositiveSmallIntegerField')()),
            ('cell_type', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'exo', ['Raw_sample'])

        # Adding model 'Gene'
        db.create_table(u'exo_gene', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('genome_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Genome'])),
            ('identifier', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('start_pos', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('end_pos', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('source', self.gf('django.db.models.fields.CharField')(max_length=100)),
        ))
        db.send_create_signal(u'exo', ['Gene'])

        # Adding unique constraint on 'Gene', fields ['genome_obj', 'identifier']
        db.create_unique(u'exo_gene', ['genome_obj_id', 'identifier'])

        # Adding unique constraint on 'Gene', fields ['genome_obj', 'name']
        db.create_unique(u'exo_gene', ['genome_obj_id', 'name'])

        # Adding model 'Peak'
        db.create_table(u'exo_peak', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('run', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Sample_run'])),
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
        ))
        db.send_create_signal(u'exo', ['Peak'])

        # Adding model 'Chromosome'
        db.create_table(u'exo_chromosome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('genome_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Genome'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('length', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'exo', ['Chromosome'])

        # Adding unique constraint on 'Chromosome', fields ['genome_obj', 'name']
        db.create_unique(u'exo_chromosome', ['genome_obj_id', 'name'])

        # Adding model 'Peak_call_parameter'
        db.create_table(u'exo_peak_call_parameter', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('run', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Sample_run'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal(u'exo', ['Peak_call_parameter'])

        # Adding unique constraint on 'Peak_call_parameter', fields ['run', 'name']
        db.create_unique(u'exo_peak_call_parameter', ['run_id', 'name'])

        # Adding model 'Sample'
        db.create_table(u'exo_sample', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('mapped_tags', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('unique_tags', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('positive_unique_count', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('negative_unique_count', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('reference_genome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Genome'])),
            ('mapping_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Mapping'])),
        ))
        db.send_create_signal(u'exo', ['Sample'])

        # Adding model 'Peak_pair_annotation'
        db.create_table(u'exo_peak_pair_annotation', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('peak_pair', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Peak_pair'])),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Gene'])),
            ('distance', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=30)),
        ))
        db.send_create_signal(u'exo', ['Peak_pair_annotation'])

        # Adding unique constraint on 'Peak_pair_annotation', fields ['peak_pair', 'gene']
        db.create_unique(u'exo_peak_pair_annotation', ['peak_pair_id', 'gene_id'])

        # Adding model 'Mapping'
        db.create_table(u'exo_mapping', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('method_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Mapping_method'])),
            ('sample_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Raw_sample'])),
        ))
        db.send_create_signal(u'exo', ['Mapping'])

        # Adding model 'Mapping_method'
        db.create_table(u'exo_mapping_method', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('mapping_method', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('mapping_method_version', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'exo', ['Mapping_method'])

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
        ))
        db.send_create_signal(u'exo', ['Peak_pair'])

        # Adding unique constraint on 'Peak_pair', fields ['pairing', 'peak1', 'peak2']
        db.create_unique(u'exo_peak_pair', ['pairing_id', 'peak1_id', 'peak2_id'])

        # Adding model 'Genome'
        db.create_table(u'exo_genome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('assembly_name', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('alternative_assembly_name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('species', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('size', self.gf('django.db.models.fields.BigIntegerField')()),
        ))
        db.send_create_signal(u'exo', ['Genome'])

        # Adding model 'Pairing_parameter'
        db.create_table(u'exo_pairing_parameter', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('pairing_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Pairing'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal(u'exo', ['Pairing_parameter'])

        # Adding unique constraint on 'Pairing_parameter', fields ['pairing_obj', 'name']
        db.create_unique(u'exo_pairing_parameter', ['pairing_obj_id', 'name'])

        # Adding model 'Pairing'
        db.create_table(u'exo_pairing', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('run', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Sample_run'])),
            ('method', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('version', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'exo', ['Pairing'])

        # Adding model 'Sample_run'
        db.create_table(u'exo_sample_run', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sample_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Sample'])),
            ('method', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('version', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'exo', ['Sample_run'])

        # Adding model 'Mapping_parameters'
        db.create_table(u'exo_mapping_parameters', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('mapping_obj', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['exo.Mapping'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal(u'exo', ['Mapping_parameters'])


    def backwards(self, orm):
        # Removing unique constraint on 'Pairing_parameter', fields ['pairing_obj', 'name']
        db.delete_unique(u'exo_pairing_parameter', ['pairing_obj_id', 'name'])

        # Removing unique constraint on 'Peak_pair', fields ['pairing', 'peak1', 'peak2']
        db.delete_unique(u'exo_peak_pair', ['pairing_id', 'peak1_id', 'peak2_id'])

        # Removing unique constraint on 'Peak_pair_annotation', fields ['peak_pair', 'gene']
        db.delete_unique(u'exo_peak_pair_annotation', ['peak_pair_id', 'gene_id'])

        # Removing unique constraint on 'Peak_call_parameter', fields ['run', 'name']
        db.delete_unique(u'exo_peak_call_parameter', ['run_id', 'name'])

        # Removing unique constraint on 'Chromosome', fields ['genome_obj', 'name']
        db.delete_unique(u'exo_chromosome', ['genome_obj_id', 'name'])

        # Removing unique constraint on 'Gene', fields ['genome_obj', 'name']
        db.delete_unique(u'exo_gene', ['genome_obj_id', 'name'])

        # Removing unique constraint on 'Gene', fields ['genome_obj', 'identifier']
        db.delete_unique(u'exo_gene', ['genome_obj_id', 'identifier'])

        # Deleting model 'Raw_sample'
        db.delete_table(u'exo_raw_sample')

        # Deleting model 'Gene'
        db.delete_table(u'exo_gene')

        # Deleting model 'Peak'
        db.delete_table(u'exo_peak')

        # Deleting model 'Chromosome'
        db.delete_table(u'exo_chromosome')

        # Deleting model 'Peak_call_parameter'
        db.delete_table(u'exo_peak_call_parameter')

        # Deleting model 'Sample'
        db.delete_table(u'exo_sample')

        # Deleting model 'Peak_pair_annotation'
        db.delete_table(u'exo_peak_pair_annotation')

        # Deleting model 'Mapping'
        db.delete_table(u'exo_mapping')

        # Deleting model 'Mapping_method'
        db.delete_table(u'exo_mapping_method')

        # Deleting model 'Peak_pair'
        db.delete_table(u'exo_peak_pair')

        # Deleting model 'Genome'
        db.delete_table(u'exo_genome')

        # Deleting model 'Pairing_parameter'
        db.delete_table(u'exo_pairing_parameter')

        # Deleting model 'Pairing'
        db.delete_table(u'exo_pairing')

        # Deleting model 'Sample_run'
        db.delete_table(u'exo_sample_run')

        # Deleting model 'Mapping_parameters'
        db.delete_table(u'exo_mapping_parameters')


    models = {
        u'exo.chromosome': {
            'Meta': {'unique_together': "(('genome_obj', 'name'),)", 'object_name': 'Chromosome'},
            'genome_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Genome']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        u'exo.gene': {
            'Meta': {'unique_together': "(('genome_obj', 'identifier'), ('genome_obj', 'name'))", 'object_name': 'Gene'},
            'end_pos': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'genome_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Genome']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'identifier': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'start_pos': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'exo.genome': {
            'Meta': {'object_name': 'Genome'},
            'alternative_assembly_name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'assembly_name': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'size': ('django.db.models.fields.BigIntegerField', [], {}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'exo.mapping': {
            'Meta': {'object_name': 'Mapping'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Mapping_method']"}),
            'sample_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Raw_sample']"})
        },
        u'exo.mapping_method': {
            'Meta': {'object_name': 'Mapping_method'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'mapping_method': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'mapping_method_version': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        u'exo.mapping_parameters': {
            'Meta': {'object_name': 'Mapping_parameters'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'mapping_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Mapping']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        },
        u'exo.pairing': {
            'Meta': {'object_name': 'Pairing'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Sample_run']"}),
            'version': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        u'exo.pairing_parameter': {
            'Meta': {'unique_together': "(('pairing_obj', 'name'),)", 'object_name': 'Pairing_parameter'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'pairing_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Pairing']"}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        },
        u'exo.peak': {
            'Meta': {'object_name': 'Peak'},
            'chrom': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'end': ('django.db.models.fields.PositiveIntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'p_value': ('django.db.models.fields.FloatField', [], {}),
            'q_value': ('django.db.models.fields.FloatField', [], {}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Sample_run']"}),
            'size': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'start': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'summit_pos_5': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'summit_pos_pileup': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'summit_val_5': ('django.db.models.fields.FloatField', [], {}),
            'summit_val_pileup': ('django.db.models.fields.FloatField', [], {}),
            'tag_count': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'exo.peak_call_parameter': {
            'Meta': {'unique_together': "(('run', 'name'),)", 'object_name': 'Peak_call_parameter'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Sample_run']"}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        },
        u'exo.peak_pair': {
            'Meta': {'unique_together': "(('pairing', 'peak1', 'peak2'),)", 'object_name': 'Peak_pair'},
            'average_height': ('django.db.models.fields.FloatField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'pairing': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Pairing']"}),
            'peak1': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['exo.Peak']"}),
            'peak2': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['exo.Peak']"}),
            'shift_distance': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'summit_height': ('django.db.models.fields.FloatField', [], {}),
            'total_score': ('django.db.models.fields.FloatField', [], {})
        },
        u'exo.peak_pair_annotation': {
            'Meta': {'unique_together': "(('peak_pair', 'gene'),)", 'object_name': 'Peak_pair_annotation'},
            'distance': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'peak_pair': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Peak_pair']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        u'exo.raw_sample': {
            'Meta': {'object_name': 'Raw_sample'},
            'antibody': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'cell_type': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'factor': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '100'}),
            'replicate_num': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'tag_length': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'total_tags': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'exo.sample': {
            'Meta': {'object_name': 'Sample'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'mapped_tags': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'mapping_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Mapping']"}),
            'negative_unique_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'positive_unique_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Genome']"}),
            'unique_tags': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'exo.sample_run': {
            'Meta': {'object_name': 'Sample_run'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'method': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'sample_obj': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Sample']"}),
            'version': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        }
    }

    complete_apps = ['exo']