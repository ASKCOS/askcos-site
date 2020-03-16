from django.urls import re_path
from rest_framework.routers import SimpleRouter

from askcos_site import api2

app_name = 'v2'

router = SimpleRouter()
router.register(r'buyables', api2.buyables.BuyablesViewSet, basename='buyables')

urlpatterns = router.urls

urlpatterns += [
    re_path(r'^context/$', api2.context.neural_network, name='context_api'),
    re_path(r'^fast-filter/$', api2.fast_filter.fast_filter, name='fast_filter_api'),
    re_path(r'^forward/$', api2.forward.template_free, name='forward_api'),
    re_path(r'^retro/$', api2.retro.singlestep, name='retro_api'),
]
