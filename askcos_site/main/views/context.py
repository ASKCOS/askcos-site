import json

from django.http import JsonResponse
from django.shortcuts import render
from django.template.loader import render_to_string

from askcos_site.askcos_celery.contextrecommender.cr_network_worker import get_n_conditions
from askcos_site.main.utils import ajax_error_wrapper


#@login_required
def context_rxnsmiles(request, smiles=None, reactants=None, product=None):
    context = {}
    if smiles is not None:
        context['reactants'] = smiles.split('>')[0]
        context['product'] = smiles.split('>')[-1]
    else:
        if reactants is not None:
            context['reactants'] = reactants
        if product is not None:
            context['product'] = product
    return render(request, 'context.html', context)

@ajax_error_wrapper
def ajax_context_rxnsmiles(request):
    '''Evaluate rxn_smiles'''
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
    verbose = json.loads(request.GET.get('verbose', 'false'))
    context_recommender = request.GET.get('context_recommender', 'Neural_Network')

    reactants = smiles.split('>>')[0].split('.')
    products = smiles.split('>>')[1].split('.')
    print('...trying to get predicted context')

    res = get_n_conditions.delay(smiles, n=10, singleSlvt=False)
    contexts = res.get(60)
    print('Got context(s)')
    print(contexts)
    if contexts is None:
        raise ValueError('Context recommender was unable to get valid context(?)')

    if contexts:
        data['html'] = render_to_string('context_recs_only.html',
            {'contexts': [context_to_dict(x) for x in contexts],
             'reactants': '.'.join(reactants)})
    else:
        data['html'] = 'No recommendations found? That is weird...'

    return JsonResponse(data)

def context_to_dict(context):
    (T1, slvt1, rgt1, cat1, t1, y1) = context
    parsible_slvt = '.'.join([x for x in slvt1.split('.') if 'Reaxys' not in x])
    nonparsible_slvt = '.'.join([x for x in slvt1.split('.') if 'Reaxys' in x])
    parsible_rgt = '.'.join([x for x in rgt1.split('.') if 'Reaxys' not in x])
    nonparsible_rgt = '.'.join([x for x in rgt1.split('.') if 'Reaxys' in x])
    parsible_cat = '.'.join([x for x in cat1.split('.') if 'Reaxys' not in x])
    nonparsible_cat = '.'.join([x for x in cat1.split('.') if 'Reaxys' in x])

    return {
        'temperature': T1,
        'solvents': parsible_slvt if parsible_slvt != '.' else '',
        'solvents_name_only': nonparsible_slvt if nonparsible_slvt != '.' else '',
        'reagents': parsible_rgt if parsible_rgt != '.' else '',
        'reagents_name_only': nonparsible_rgt if nonparsible_rgt != '.' else '',
        'reagents_combined': '.'.join([x for x in parsible_rgt.split('.') + parsible_cat.split('.') if x]),
        'catalysts': parsible_cat if parsible_cat != '.' else '',
        'catalysts_name_only': nonparsible_cat if nonparsible_cat != '.' else '',
        'time': t1,
        'yield': y1,
    }

#@login_required
def context_rxnsmiles_target(request, smiles):
    '''Synth interactive initialized w/ reaction smiles'''
    return context_rxnsmiles(request, reactants=smiles.split('>')[0], product=smiles.split('>')[-1])

#@login_required
def context_rxnsmiles_target2(request, reactants, product):
    '''Synth interactive initialized w/ reaction smiles'''
    return context_rxnsmiles(request, reactants=reactants, product=product)
