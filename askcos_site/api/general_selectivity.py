from rdkit import Chem
from django.http import JsonResponse
from celery.exceptions import TimeoutError
from askcos_site.askcos_celery.generalselectivity.selec_worker import get_selec

TIMEOUT = 60

def selectivity(request):
    '''Evaluate rxn_smiles'''
    resp = {}
    resp['request'] = dict(**request.GET)
    run_async = request.GET.get('async', False)
    rxn_smiles = request.GET.get('rxn_smiles', None)

    if rxn_smiles is None:
        resp['error'] = 'Required target "rxn_smiles" missing'
        return JsonResponse(resp, status=400)

    r, _, p = rxn_smiles.split('>')
    r = Chem.MolFromSmiles(r)
    p = Chem.MolFromSmiles(p)

    if not r or not p:
        resp['error'] = 'Cannot parse smiles with rdkit'
        return JsonResponse(resp, status=400)

    res = get_selec.delay(rxn_smiles)
    try:
        result = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request time out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)

    resp['results'] = {'scores': result}

    return JsonResponse(resp)

