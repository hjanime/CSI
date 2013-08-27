from django.contrib import admin
from exo.models import *

admin.site.register(Genome)
#admin.site.register(Raw_sample)
admin.site.register(Sample)

admin.site.register(Peak)

admin.site.register(Method)
admin.site.register(Mapping)

admin.site.register(Call_peak)
admin.site.register(Pairing)

admin.site.register(Chromosome)
admin.site.register(Gene)


class ParameterInline(admin.StackedInline):
    model = Parameter
    extra = 0
    


class MappingInline(admin.StackedInline):
    model = Mapping
    extra = 0
    
class RSAdmin(admin.ModelAdmin):
    inlines = [MappingInline]
    extra = 0
    
class PSAdmin(admin.ModelAdmin):
    inlines = [ParameterInline, MappingInline]

admin.site.register(Raw_sample, RSAdmin )

#admin.site.register(Parameter)
admin.site.register(Parameter_set, PSAdmin)

