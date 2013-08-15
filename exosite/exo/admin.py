from django.contrib import admin
from exo.models import *

admin.site.register(Genome)
admin.site.register(Raw_sample)
admin.site.register(Sample)
admin.site.register(Parameter_set)
admin.site.register(Peak)

admin.site.register(Method)
admin.site.register(Mapping)
admin.site.register(Parameter)
admin.site.register(Call_peak)
admin.site.register(Pairing)

admin.site.register(Chromosome)
admin.site.register(Gene)
