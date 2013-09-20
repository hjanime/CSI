'''
Created on 28 Aug, 2013

@author: caofan
'''
from django import template
register = template.Library()


@register.inclusion_tag('peak_rows.html')
def load_peaks(peaks, page):
    return {'peak_list':peaks, "page": page}