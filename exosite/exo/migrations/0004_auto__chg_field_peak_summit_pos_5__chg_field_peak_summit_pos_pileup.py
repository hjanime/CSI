# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):

        # Changing field 'Peak.summit_pos_5'
        db.alter_column(u'exo_peak', 'summit_pos_5', self.gf('django.db.models.fields.IntegerField')())

        # Changing field 'Peak.summit_pos_pileup'
        db.alter_column(u'exo_peak', 'summit_pos_pileup', self.gf('django.db.models.fields.IntegerField')())

    def backwards(self, orm):

        # Changing field 'Peak.summit_pos_5'
        db.alter_column(u'exo_peak', 'summit_pos_5', self.gf('django.db.models.fields.PositiveIntegerField')())

        # Changing field 'Peak.summit_pos_pileup'
        db.alter_column(u'exo_peak', 'summit_pos_pileup', self.gf('django.db.models.fields.PositiveIntegerField')())

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
            'fold': ('django.db.models.fields.FloatField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'p_value': ('django.db.models.fields.FloatField', [], {}),
            'q_value': ('django.db.models.fields.FloatField', [], {}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['exo.Call_peak']"}),
            'size': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'start': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'summit_pos_5': ('django.db.models.fields.IntegerField', [], {}),
            'summit_pos_pileup': ('django.db.models.fields.IntegerField', [], {}),
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
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '100'}),
            'negative_unique_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'positive_unique_count': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'unique_tags': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'updated_at': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['exo']