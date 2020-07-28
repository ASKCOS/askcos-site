from django.shortcuts import render


def synth_interactive(request):
    context = {}
    return render(request, 'forward_interactive.html', context)
