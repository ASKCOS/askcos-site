from bson.objectid import ObjectId
from django.shortcuts import render, HttpResponse

from askcos_site.globals import retro_templates


def template_target_export(request, id):
    rx_ids = template_target(request, id, return_refs_only=True)

    txt = '{"fileName":"reaxys_query.json", "version":"1.0","content": {"id":"root","facts": ['
    for i, rx_id in enumerate(rx_ids):
        if i != 0:
            txt += ','
        txt += '{"id":"Reaxys487",'
        if i != 0:
            txt += '"logicOperator":"OR",'
        txt += '"fields": [{"value":"%i","boundOperator":"op_num_equal","id":"RX.ID","displayName":"Reaction ID"}],' % (rx_id)
        txt += '"fieldsLogicOperator":"AND","exist":false,"bio":false}'

    txt += '], "exist":false,"bio":false } }'

    response = HttpResponse(txt, content_type='text/csv')
    response['Content-Disposition'] = 'attachment;filename=reaxys_query.json'
    return response

def template_view(request):
    return render(request, 'template_view.html')

#@login_required
def template_target(request, id, return_refs_only=False):
    '''
    Examines a template from its id in the database
    where id == str(_id)

    Reaxys templates refer to instances, but we need the
    reactions for visualization

    If return_refs_only, builds .json query for Reaxys
    '''
    context = {}

    transform = retro_templates.find_one({'_id': id})
    if not transform:
        try:
            transform = retro_templates.find_one({'_id': ObjectId(id)})
        except:
            transform = None

    if not transform:
        context['err'] = 'Transform not found'
        return render(request, 'template.html', context)

    context['template'] = transform
    context['template']['id'] = id
    reference_ids = transform['references']

    if return_refs_only:
        return list(sorted(set(int(_id.split('-')[0]) for _id in reference_ids)))[:500]

    references = []
    rx_docs = {}; xrn_to_smiles = {}
    context['total_references'] = len(reference_ids)

    context['cannot_view_any'] = True
    context['ref_ids'] = '; '.join([str(y) for y in sorted(set(int(_id.split('-')[0]) for _id in reference_ids))])

    return render(request, 'template.html', context)
