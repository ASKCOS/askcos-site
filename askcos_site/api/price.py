from django.http import JsonResponse

from askcos_site.globals import pricer

def price(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    smiles = request.GET.get('smiles', None)
    try:
        p = pricer.lookup_smiles(smiles, alreadyCanonical=True)
    except Exception as e:
        resp['error'] = str(e)
        return JsonResponse(resp, status=400)
    resp['price'] = p
    return JsonResponse(resp)
