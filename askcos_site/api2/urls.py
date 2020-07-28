from django.urls import path
from rest_framework.routers import SimpleRouter
from rest_framework_jwt.views import obtain_jwt_token, refresh_jwt_token

from askcos_site import api2

app_name = 'v2'

router = SimpleRouter()
router.register(r'buyables', api2.buyables.BuyablesViewSet, basename='buyables_api')
router.register(r'rdkit/smiles', api2.rdkit.SmilesViewSet, basename='smiles_api')
router.register(r'template', api2.template.TemplateViewSet, basename='template_api')
router.register(r'results', api2.results.ResultsViewSet, basename='results_api')
router.register(r'celery/task', api2.celery.CeleryTaskViewSet, basename='celery_task_api')
router.register(r'banlist/chemicals', api2.banlist.BannedChemicalsViewSet, basename='banlist_chemicals_api')
router.register(r'banlist/reactions', api2.banlist.BannedReactionsViewSet, basename='banlist_reactions_api')

urlpatterns = router.urls

urlpatterns += [
    path('atom-mapper/', api2.atom_mapper.atom_mapper, name='atom_mapper_api'),
    path('celery/', api2.celery.celery_status, name='celery_api'),
    path('cluster/', api2.cluster.cluster, name='cluster_api'),
    path('context/', api2.context.neural_network, name='context_api'),
    path('draw/', api2.draw.drawer, name='draw_api'),
    path('fast-filter/', api2.fast_filter.fast_filter, name='fast_filter_api'),
    path('forward/', api2.forward.template_free, name='forward_api'),
    path('impurity/', api2.impurity.impurity_predict, name='impurity_api'),
    path('reactions/', api2.reactions.reactions, name='reactions_api'),
    path('retro/', api2.retro.singlestep, name='retro_api'),
    path('retro/models/', api2.retro.models, name='retro_models_api'),
    path('scscore/', api2.scscore.scscore, name='scscore_api'),
    path('selectivity/', api2.selectivity.selectivity, name='selectivity_api'),
    path('general-selectivity/', api2.general_selectivity.selectivity, name='general_selectivity_api'),
    path('tree-builder/', api2.tree_builder.tree_builder, name='tree_builder_api'),

    path('token-auth/', obtain_jwt_token, name='token_auth_api'),
    path('token-refresh/', refresh_jwt_token, name='token_refresh_api'),
]

urlpatterns.append(path('', api2.root.RootAPIView.as_view(namespace=app_name, urlpatterns=urlpatterns), name='root_api'))
