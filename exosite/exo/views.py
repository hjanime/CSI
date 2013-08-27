# Create your views here.
from django.http import HttpResponse
from django.shortcuts import render, render_to_response
from django.template import Context

from exo.models import *

def index(request):
    return render(request, 'index.html')

def raw(request):
    raw_samples = Raw_sample.objects.order_by('name')
    headers = ["Sample name", "Cell line", "Factor",
               "Antibody", "Tag length", "Tag count",
               "Coverage", "Replicate#", "Linked mappings"]
    item_list = []
    mapping_sets = []
    curr_ids = []
    for r in raw_samples:
        temp_item = [r.name, r.cell_type, r.factor, 
                     r.antibody, r.tag_length,
                     r.total_tags, r.coverage, r.replicate_num]
        item_list.append(temp_item)
        mapping_sets.append(r.mapping_set.all())
        curr_ids.append( r.id )
        
    context = Context({"item_list": zip(item_list, mapping_sets, curr_ids),
                       "table_header": headers })
    return render(request, 'raw.html', context)

def mapped_helper(request, samples):
    '''
    This is a helper function for processing mapped view
    '''
    headers = ["Name", "Mapped tags", "Unique tags", "Positive unique tags", "Negative unique tags",
                "Raw sample", "Genome assembly", "Method", "Peak sets"]
    item_list = []
    raw_samples = []
    genomes = []
    params = []
    peak_sets = []
    for m in samples:
        mapping_obj = m.mapping_obj
        temp_peak_sets = m.call_peak_set.all()
        temp_tuple = ([m.name, m.mapped_tags, m.unique_tags, 
                       m.positive_unique_count, m.negative_unique_count],
                       mapping_obj.sample_obj, mapping_obj.genome_obj, mapping_obj.params,temp_peak_sets)
        item_list.append( temp_tuple )
    context = Context({"item_list": item_list,
                       "table_header": headers})
    
    return render(request, 'map.html', context)

def mapped(request):
    samples = Sample.objects.order_by("name")
    return mapped_helper(request, samples)
        
def mapped_for_raw(request, raw_id):
    mappings = Mapping.objects.filter(sample_obj__in=[raw_id,])
    samples = Sample.objects.filter(mapping_obj__in=mappings)
    return mapped_helper(request, samples)

def peaks(request):
    return HttpResponse("Hello")


def peakpairs(request):
    return HttpResponse("Hello")


def genomes(request):
    return HttpResponse("Hello")


def chroms(request):
    return HttpResponse("Hello")

def genes(request):
    return HttpResponse("Hello")
