import time

from django.http import JsonResponse
from django.shortcuts import render
from django.urls import reverse

from askcos_site.askcos_celery.atom_mapper.atom_mapping_worker import get_atom_mapping
from ..utils import ajax_error_wrapper


# @login_required
def atom_mapping(request):
    return render(request, 'mapping.html')

@ajax_error_wrapper
def ajax_find_atom_mapping(request):
    '''Perform the forward synthesis'''
    data = {'err': False}

    rxnsmiles = request.GET.get('rxnsmiles', '')

    mapper = request.GET.get('mapper', 'WLN atom mapper')

    print('reactants: {}'.format(rxnsmiles))
    print('mapper: {}'.format(mapper))

    startTime = time.time()

    ## need to work on it
    result = get_atom_mapping.delay(rxnsmiles, mapper=mapper)
    data['rxnsmiles_mapped'] = result.get(10)
    print('---------------------------------')
    print(data)
    print('---------------------------------')

    data['html_time'] = '{:.3f} seconds elapsed'.format(time.time() - startTime)

    if data['rxnsmiles_mapped']:
        url_mapping = reverse('draw_mapped_reaction', kwargs={'smiles': data['rxnsmiles_mapped']})
        url_highlight = reverse('draw_highlighted_reaction', kwargs={'smiles': data['rxnsmiles_mapped']})
        data['html_mapping'] = '<img src="' + url_mapping + '">'
        data['html_highlight'] = '<img src="' + url_highlight + '">'
    else:
        data['html_mapping'] = 'No atom mapping is found'
        data['html_highlight'] = 'No highlighted reaction iamge to show'
    return JsonResponse(data)
