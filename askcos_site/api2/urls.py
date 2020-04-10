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
router.register(r'blacklist/chemicals', api2.blacklist.BlacklistedChemicalsViewSet, basename='blacklist_chemicals_api')
router.register(r'blacklist/reactions', api2.blacklist.BlacklistedReactionsViewSet, basename='blacklist_reactions_api')

urlpatterns = router.urls

urlpatterns += [
    path('atom-mapper/', api2.atom_mapper.atom_mapper, name='atom_mapper_api'),
    path('celery/', api2.celery.celery_status, name='celery_api'),
    path('cluster/', api2.cluster.cluster, name='cluster_api'),
    path('context/', api2.context.neural_network, name='context_api'),
    path('fast-filter/', api2.fast_filter.fast_filter, name='fast_filter_api'),
    path('forward/', api2.forward.template_free, name='forward_api'),
    path('impurity/', api2.impurity.impurity_predict, name='impurity_api'),
    path('reactions/', api2.reactions.reactions, name='reactions_api'),
    path('retro/', api2.retro.singlestep, name='retro_api'),
    path('scscore/', api2.scscore.scscore, name='scscore_api'),
    path('selectivity/', api2.selectivity.selectivity, name='selectivity_api'),
    path('tree-builder/', api2.tree_builder.tree_builder, name='tree_builder_api'),

    path('token-auth/', obtain_jwt_token, name='token_auth_api'),
    path('token-refresh/', refresh_jwt_token, name='token_refresh_api'),

    path('', api2.root.apiroot, name='root_api'),
]
