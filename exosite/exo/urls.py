'''
Created on 21 Aug, 2013

@author: caofan
'''

from django.conf.urls import patterns, url
from exo import views

urlpatterns = patterns('', 
                       url(r'^$', views.index, name='index'),
                       url(r'^raw/$', views.raw, name='raw'),
                       url(r'^mapped/$', views.mapped, name='mapped'),
                       
                       url(r'^peaks/$', views.peaks, name="peaks"),
                       url(r'^peakpairs/$', views.peakpairs, name="peakpairs"),
                       url(r'^genomes/$', views.genomes, name="genomes"),
                       url(r'^chroms/$', views.chroms, name="chroms"),
                       url(r'^genes/$', views.genes, name="genes"),
                       url(r'^mapped/for_raw/(?P<raw_id>\d+)/$',views.mapped_for_raw, name="mapped_for_raw"),
                       url(r'^peaks/for_call/(?P<call_id>\d+)/$', views.peak_view, name="peaks_for_call"),
                       )