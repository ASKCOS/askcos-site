import json

from django.http import JsonResponse
from django.shortcuts import render, HttpResponse

from askcos_site.globals import pricer
from askcos_site.main.views.users import can_modify_buyables
from askcos_site.main.utils import ajax_error_wrapper


def price_smiles_func(smiles):
    return pricer.lookup_smiles(smiles, alreadyCanonical=True)

def price_smiles(request, smiles):
    response = HttpResponse(content_type = 'text/plain')
    ppg = pricer.lookup_smiles(smiles, alreadyCanonical=True)
    if ppg:
        response.write('${}/g'.format(ppg))
    else:
        response.write('cannot buy')
    return response

@ajax_error_wrapper
def ajax_price_smiles(request):
    print('Got price request')
    smiles = request.GET.get('smiles', None)
    data = {'err': False, 'smiles': smiles}
    isomericSmiles = json.loads(request.GET.get('isomericSmiles', 'false'))
    print('isomericSmiles: {}'.format(isomericSmiles))
    data['ppg'] = pricer.lookup_smiles(smiles, alreadyCanonical=False, isomericSmiles=isomericSmiles)
    print('Result: {}'.format(data['ppg']))
    data['buyable'] = data['ppg'] != 0.
    if data['ppg'] == 0.:
        data['html'] = 'This chemical is <b>not</b> in our database currently'
    else:
        data['html'] = 'This chemical is purchaseable for an estimated <b>$%i/g</b>' % data['ppg']
    return JsonResponse(data)


def pricing(request):
    return render(request, 'pricing.html', {})


def buyables(request):
    return render(request, 'buyables.html', {'can_modify_buyables': can_modify_buyables(request)})


def price_xrn(request, xrn):
    response = HttpResponse(content_type = 'text/plain')
    ppg = pricer.lookup_xrn(xrn)
    if ppg:
        response.write('${}/g'.format(ppg))
    else:
        response.write('cannot buy')
    return response
