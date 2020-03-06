from django.http import JsonResponse
from rdkit import Chem

from askcos_site.globals import scscorer


def scscore(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    smiles = request.GET.get('smiles', None)

    if not smiles:
        resp['error'] = 'Required parameter "smiles" missing'
        return JsonResponse(resp, status=400)

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        resp['error'] = 'Cannot parse smiles smiles with rdkit'
        return JsonResponse(resp, status=400)

    scscore = scscorer.get_score_from_smiles(smiles, noprice=True)
    resp['score'] = scscore
    return JsonResponse(resp)
