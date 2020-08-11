from django.shortcuts import render

def ketcher_iframe(request):
    return render(request, 'ketcher_iframe.html')

def ketcher_iframe_min(request):
    return render(request, 'ketcher_iframe_min.html')
