# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting field 'Genome.id'
        db.delete_column(u'exo_genome', u'id')


        # Changing field 'Genome.assembly_name'
        db.alter_column(u'exo_genome', 'assembly_name', self.gf('django.db.models.fields.CharField')(max_length=20, primary_key=True))
        # Adding unique constraint on 'Genome', fields ['assembly_name']
        db.create_unique(u'exo_genome', ['assembly_name'])


    def backwards(self, orm):
        # Removing unique constraint on 'Genome', fields ['assembly_name']
        db.delete_unique(u'exo_genome', ['assembly_name'])


        # User chose to not deal with backwards NULL issues for 'Genome.id'
        raise RuntimeError("Cannot reverse this migration. 'Genome.id' and its values cannot be restored.")

        # Changing field 'Genome.assembly_name'
        db.alter_column(u'exo_genome', 'assembly_name', self.gf('django.db.models.fields.CharField')(max_length=20))

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
            'assembly_name': ('django.db.models.fields.CharField', [], {'max_length': '20', 'primary_key': 'True'}),
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