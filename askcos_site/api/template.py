from bson.objectid import ObjectId
from django.http import JsonResponse

from askcos_site.globals import retro_templates


def template(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    _id = request.GET.get('id')
    transform = retro_templates.find_one({'_id': _id})
    if not transform:
        transform = retro_templates.find_one({'_id': ObjectId(_id)})
    if not transform:
        resp['error'] = 'Cannot find template with id '+_id
        return JsonResponse(resp)
    transform['_id'] = _id
    transform.pop('product_smiles', None)
    transform.pop('name', None)
    refs = transform.pop('references', [''])
    transform['references'] = list(map(lambda x: x.split('-')[0], refs))
    resp['template'] = transform
    return JsonResponse(resp)
