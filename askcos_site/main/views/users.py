import os
from datetime import datetime

import django.contrib.auth.views
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.shortcuts import render

from askcos_site.globals import db_client
from askcos_site.main.models import SavedResults

results_collection = db_client['results']['results']

AUTH_MODIFY_BUYABLES = os.environ.get('AUTH_MODIFY_BUYABLES') == 'True'

can_control_robot = lambda request: request.user.get_username() in ['ccoley']

def can_view_reaxys(request):
    return request.user.is_authenticated and request.user.groups.filter(name='reaxys_view').exists()

def can_avoid_banned_chemicals(request):
    return request.user.is_authenticated and request.user.groups.filter(name='avoid_banned_chemicals').exists()

def can_modify_buyables(request):
    if not AUTH_MODIFY_BUYABLES:
        return True
    return request.user.is_authenticated and request.user.groups.filter(name='modify_buyables').exists()

def log_this_request(method):
    def f(*args, **kwargs):
        try:
            print('User %s requested view %s with args %r and kwargs %r' % \
                args[0].user.get_username(), method.__name__, args, kwargs)
        except Exception as e:
            print(e)
        return method(*args, **kwargs)
    return f


def login(request):
    '''
    User login
    '''
    return django.contrib.auth.views.login(request, template_name='login.html')

def logout(request):
    '''
    User logout
    '''
    return django.contrib.auth.views.logout(request, template_name='logout.html')

@login_required
def user_saved_results(request, err=None):
    saved_results = SavedResults.objects.filter(user=request.user)
    return render(request, 'saved_results.html', {'saved_results':saved_results, 'err': err})

@login_required
def user_saved_results_id(request, _id=-1):
    result = results_collection.find_one({'_id': _id})
    if not result:
        return user_saved_results(request, err='Could not find that ID')
    return render(request, 'saved_results_id.html',
        {'html': result['result']})

@login_required
def user_saved_results_del(request, _id=-1):
    obj = SavedResults.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        os.remove(obj[0].fpath)
        obj[0].delete()
    return user_saved_results(request, err=None)

@login_required
def ajax_user_save_page(request):
    html = request.POST.get('html', None)
    desc = request.POST.get('desc', 'no description')
    dt   = request.POST.get('datetime', datetime.utcnow().strftime('%B %d, %Y %H:%M:%S %p UTC'))
    if html is None:
        print('Got None html')
        data = {'err': 'Could not get HTML to save'}
        return JsonResponse(data)
    print('Got request to save a page')
    now = datetime.now()
    result_id = str(hash((now, request.user)))
    obj = SavedResults.objects.create(user=request.user,
        description=desc,
        dt=dt,
        created=now,
        result_id=result_id,
        result_state='N/A',
        result_type='html'
    )
    results_collection.insert_one({
        '_id': result_id,
        'result': html,
    })
    print('Created saved object {}'.format(obj.id))
    return JsonResponse({'err': False})


@login_required
def banlist(request):
    return render(request, 'banlist.html')
