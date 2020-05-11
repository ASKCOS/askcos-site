from django.http import JsonResponse
from rdkit import Chem

from askcos_site.askcos_celery.impurity.impurity_worker import get_impurities

TIMEOUT = 30


def impurity_predict(request):
    resp = {}
    resp['request'] = dict(**request.GET)

    reactants = request.GET.get('reactants')
    solvents = request.GET.get('solvents', '')
    reagents = request.GET.get('reagents', '')
    products = request.GET.get('products', '')
    top_k = int(request.GET.get('top_k', 3))
    threshold = float(request.GET.get('threshold', 0.75))
    predictor = request.GET.get('predictor', 'WLN forward predictor')
    inspector = request.GET.get('inspector', 'Reaxys inspector')
    mapper = request.GET.get('mapper', 'WLN atom mapper')
    check_mapping = request.GET.get('check_mapping', 'True') in ['True', 'true']

    # TODO: add model selections

    if not reactants:
        resp['error'] = 'Required parameter "reactants" missing'
        return JsonResponse(resp, status=400)

    rmol = Chem.MolFromSmiles(reactants)
    if not rmol:
        resp['error'] = 'Cannot parse reactants smiles with rdkit'
        return JsonResponse(resp, status=400)

    smol = Chem.MolFromSmiles(solvents)
    if not smol:
        resp['error'] = 'Cannot parse solvent smiles with rdkit'
        return JsonResponse(resp, status=400)

    remol = Chem.MolFromSmiles(reagents)
    if not remol:
        resp['error'] = 'Cannot parse reagents smiles with rdkit'
        return JsonResponse(resp, status=400)

    pmol = Chem.MolFromSmiles(products)
    if not pmol:
        resp['error'] = 'Cannot parse products smiles with rdkit'
        return JsonResponse(resp, status=400)

    # TODO: need work
    res = get_impurities.delay(
        reactants, reagents=reagents, products=products, solvents=solvents,
        top_k=top_k, threshold=threshold,
        predictor_selection=predictor,
        inspector_selection=inspector,
        mapper_selection=mapper,
        check_mapping=check_mapping
    )

    resp['task_id'] = res.task_id
    return JsonResponse(resp)
