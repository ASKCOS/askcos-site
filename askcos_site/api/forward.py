from celery.exceptions import TimeoutError
from django.http import JsonResponse
from rdkit import Chem

from askcos_site.askcos_celery.treeevaluator.template_free_forward_predictor_worker import get_outcomes

TIMEOUT = 30


def template_free(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    reactants = request.GET.get('reactants')
    solvent = request.GET.get('solvent', '')
    reagents = request.GET.get('reagents', '')
    num_results = int(request.GET.get('num_results', 100))

    if not reactants:
        resp['error'] = 'Required parameter "reactants" missing'
        return JsonResponse(resp, status=400)

    rmol = Chem.MolFromSmiles(reactants)
    if not rmol:
        resp['error'] = 'Cannot parse reactants smiles with rdkit'
        return JsonResponse(resp, status=400)

    if solvent:
        smol = Chem.MolFromSmiles(solvent)
        if not smol:
            resp['error'] = 'Cannot parse solvent smiles with rdkit'
            return JsonResponse(resp, status=400)

    if reagents:
        remol = Chem.MolFromSmiles(reagents)
        if not remol:
            resp['error'] = 'Cannot parse reagents smiles with rdkit'
            return JsonResponse(resp, status=400)

    combined_smiles = Chem.MolToSmiles(rmol, isomericSmiles=True)
    if reagents:
        combined_smiles += '.{}'.format(Chem.MolToSmiles(remol, isomericSmiles=True))
    if solvent:
        combined_smiles += '.{}'.format(Chem.MolToSmiles(smol, isomericSmiles=True))

    res = get_outcomes.delay(combined_smiles, top_n=num_results)
    
    try:
        outcomes = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)

    resp['outcomes'] = outcomes
    return JsonResponse(resp)
