from django.urls import re_path
from rest_framework.routers import SimpleRouter

from askcos_site import api2

app_name = 'v2'

router = SimpleRouter()
router.register(r'buyables', api2.buyables.BuyablesViewSet, basename='buyables')

urlpatterns = router.urls

urlpatterns += [
    re_path(r'^retro/$', api2.retro.singlestep, name='retro_api'),
]
