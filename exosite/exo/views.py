# Create your views here.
from django.http import HttpResponse, Http404
from django.shortcuts import render, render_to_response
from django.template import Context, RequestContext, loader
from django.core.paginator import Paginator, InvalidPage
import sys

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
    if request.GET.get('map_id'):
        samples = Sample.objects.filter(id=request.GET.get('map_id'))
    return mapped_helper(request, samples)
        
def mapped_for_raw(request, raw_id):
    mappings = Mapping.objects.filter(sample_obj__in=[raw_id,])
    samples = Sample.objects.filter(mapping_obj__in=mappings)
    return mapped_helper(request, samples)

def peaks(request):
    '''
    List out all peak sets.
    '''
    calls = Call_peak.objects.all()
    headers = ["Sample", "# Peaks", "Method", '']
    item_list = []
    for c in calls:
        #item_list.append((c.sample_obj, c.peak_set.count(), c.params))
        item_list.append(c)
    context= Context({"item_list": item_list,
                      "table_header": headers})
    return render(request, 'peaksets.html', context)

def peaks_1_index(request, call_id):
    '''
    Pull the first page of all peaks of one set.
    '''
    run_obj = Call_peak.objects.get(id=call_id)
    peaks = Peak.objects.filter(run=call_id)
    paginator = Paginator(peaks, 50)
    try:
        page_obj = paginator.page(1)
    except InvalidPage:
        raise Http404
    headers = ["Chrom","Start","End","Size","Strand","Tag count","Summit positionP",
               "Summit valueP", "Summit pos5", "Summit value5", "P-value", "q-value", 'Fold']
    context = Context({"peak_list": page_obj.object_list,
                       "page": page_obj,
                       "table_header": headers,
                       "run_obj": run_obj})
    
    return render(request, "peaks.html", context)

def peaks_paginator_ajax(request, call_id, page_num):
    '''
    Ajax request for more data upon scrolling.
    '''
    peaks = Peak.objects.filter(run=call_id)
    
    paginator = Paginator(peaks, 50)
    
    try:
        page_obj = paginator.page(page_num)
    except InvalidPage:
        raise Http404
    context = Context({"peak_list": page_obj.object_list,
                       "page": page_obj,
                       })
    return render(request, "peak_rows.html",context)


def __get_float_range(request, minvname, maxvname, vrange):
    if request.POST.get(minvname):
        try:
            minv = float(request.POST.get(minvname))
            vrange[0] = minv
        except ValueError:
            pass
        try:
            maxv = float(request.POST.get(maxvname))
            vrange[1] = maxv
        except ValueError:
            pass
    return vrange
            

def __retrieve_and_apply_filter_from_post(request, querySet):
    filter_size = [-100, sys.maxint]
    filter_pvalue = [0, sys.maxint]
    filter_qvalue = [0, sys.maxint]
    filter_fold = [0, sys.maxint]
    filter_tagcount = [0, sys.maxint]
    filter_summitvp = [0, sys.maxint]
    filter_summitv5 = [0, sys.maxint]
    
    filter_size = __get_float_range(request,'filter_minsize', 'filter_maxsize', filter_size)
    filter_pvalue = __get_float_range(request,'filter_minpvalue', 'filter_maxpvalue', filter_pvalue)
    filter_qvalue = __get_float_range(request,'filter_minqvalue', 'filter_maxqvalue', filter_qvalue)
    filter_fold = __get_float_range(request,'filter_minfold', 'filter_maxfold', filter_fold)
    filter_summitvp = __get_float_range(request,'filter_minsummitvp', 'filter_maxsummitvp', filter_summitvp)
    filter_summitv5 = __get_float_range(request,'filter_minsummitv5', 'filter_maxsummitv5', filter_summitv5)
    filter_tagcount = __get_float_range(request, 'filter_mintagcount', 'filter_maxtagcount', filter_tagcount)
    if request.POST.get('filter_chrom'):
        filter_chrom = set([ i.strip() for i in request.POST.get('filter_chrom').split(',') ])
        if '' in filter_chrom:
            filter_chrom.remove('')
        filter_chrom = list(filter_chrom)
        if len(filter_chrom) > 0:
            querySet = querySet.filter(chrom__in=filter_chrom)
    return querySet.filter(size__range = filter_size,
                               p_value__range = filter_pvalue,
                               q_value__range = filter_qvalue,
                               fold__range = filter_fold,
                               summit_val_pileup__range = filter_summitvp,
                               summit_val_5__range = filter_summitv5,
                               tag_count__range = filter_tagcount)



def __retrieve_and_apply_sort_from_post(request, querySet):
    if request.POST.get('sort_attr') and request.POST.get('sort_direction'):
        sort_attrs = [i.strip().lower() for i in request.POST.getlist('sort_attr')]
        sort_directs = [i.strip().lower() for i in request.POST.getlist('sort_direction')]
        
        ATTR_MAP_TO_DB = {'location': ['chrom','start','end'],
                          'strand': ['strand',],
                          'pvalue': ['p_value'],
                          'qvalue': ['q_value'],
                          'fold': ['fold'],
                          'tagcount': ['tag_count'],
                          'summitvp': ['summit_val_pileup'],
                          'summitv5': ['summit_val_5'],
                          'size': ['size'],
                          }
        db_attrs = []
        for i, attr in enumerate(sort_attrs):
            temp = ATTR_MAP_TO_DB[attr]
            if i < len(sort_directs) and sort_directs[i] == 'desc':
                temp = ['-' + k for k in temp]
            db_attrs += temp
        if len(db_attrs) > 0:
            querySet = querySet.order_by(*db_attrs)    
    return querySet

def peak_view(request, call_id):
    
    headers = ["Chrom","Start","End","Size","Strand","Tag count","Summit positionP",
               "Summit valueP", "Summit pos5", "Summit value5", "P-value", "q-value", 'Fold']
    context = Context({"table_header": headers,
                       "run_obj": Call_peak.objects.get(id=call_id)})
    
    peaks = Peak.objects.filter(run=call_id)
    page_num = 1
    template = 'peaks.html'
    if request.method == 'GET':
        '''
        If request uses GET method, assume no filter
        and sorting need to be applied.
        '''
        paginator = Paginator(peaks, 50)
        
        if request.GET.get('page_num'):
            page_num = request.GET.get('page_num')
            template = 'peak_rows.html'
        
        context['has_filter'] = False
    elif request.method == 'POST':
        '''
        Need to retrieve all data from request
        by POST method.
        '''
        peaks = __retrieve_and_apply_filter_from_post(request, peaks)
        peaks = __retrieve_and_apply_sort_from_post(request, peaks)
        paginator = Paginator(peaks, 50)
        
        if request.POST.get('page_num'):
            page_num = request.POST.get('page_num')
            template = 'peak_rows.html'
        context['has_filter'] = True
    try:
        page_obj = paginator.page(page_num)
    except InvalidPage:
        raise Http404
    
    context["peak_list"] = page_obj.object_list
    context["page"] = page_obj
    
    return render(request, template, context)
        

def peakpairs(request):
    return HttpResponse("Hello")


def genomes(request):
    return HttpResponse("Hello")


def chroms(request):
    return HttpResponse("Hello")

def genes(request):
    return HttpResponse("Hello")
