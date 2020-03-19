from django.urls import re_path
from rest_framework.routers import SimpleRouter
from rest_framework_jwt.views import obtain_jwt_token

from askcos_site import api2

app_name = 'v2'

router = SimpleRouter()
router.register(r'buyables', api2.buyables.BuyablesViewSet, basename='buyables')
router.register(r'rdkit/smiles', api2.rdkit.SmilesViewSet, basename='smiles')
router.register(r'template', api2.template.TemplateViewSet, basename='template')
router.register(r'results', api2.results.ResultsViewSet, basename='results')

urlpatterns = router.urls

urlpatterns += [
    re_path(r'^celery/$', api2.status.celery_status, name='celery_api'),
    re_path(r'^celery/task/$', api2.status.task_status, name='celery_task_api'),
    re_path(r'^cluster/$', api2.cluster.cluster, name='cluster_api'),
    re_path(r'^context/$', api2.context.neural_network, name='context_api'),
    re_path(r'^fast-filter/$', api2.fast_filter.fast_filter, name='fast_filter_api'),
    re_path(r'^forward/$', api2.forward.template_free, name='forward_api'),
    re_path(r'^impurity/$', api2.impurity.impurity_predict, name='impurity_api'),
    re_path(r'^reactions/$', api2.reactions.reactions, name='reactions_api'),
    re_path(r'^retro/$', api2.retro.singlestep, name='retro_api'),
    re_path(r'^scscore/$', api2.scscore.scscore, name='scscore_api'),
    re_path(r'^selectivity/$', api2.selectivity.selectivity, name='selectivity'),
    re_path(r'^treebuilder/$', api2.tree_builder.tree_builder, name='tree_builder_api'),

    re_path(r'^token-auth/$', obtain_jwt_token, name='token_api'),
]
