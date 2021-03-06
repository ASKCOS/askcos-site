from django.contrib.auth.decorators import login_required
from django.shortcuts import render

from askcos_site.globals import db_client
from askcos_site.main.models import SavedResults

results_collection = db_client['results']['results']


@login_required
def my_results(request):
    context = {}
    results = SavedResults.objects.filter(user=request.user)
    if results:
        context['results'] = results
    return render(request, 'results.html', context)

@login_required
def view_result(request):
    id_ = request.GET.get('id')
    try:
        result = SavedResults.objects.get(result_id=id_, user=request.user)
    except:
        return render(
            request,
            'trees_only.html',
            {'trees': []}
        )
    if result and result.result_state == 'completed':
        result_doc = results_collection.find_one({'_id': id_})
        (tree_status, trees) = result_doc['result']
        (num_chemicals, num_reactions, _) = tree_status
        return render(
            request,
            'trees_only.html',
            {'trees': trees}
        )
    else:
        return render(
            request,
            'trees_only.html',
            {'trees': []}
        )

@login_required
def view_tree_graph(request):
    id_ = request.GET.get('id')
    return render(request, 'trees_graph.html', {'id': id_})
