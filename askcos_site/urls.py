from django.urls import include, path, re_path
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
admin.autodiscover()
import django.contrib.auth.urls
from django.views.generic import TemplateView
import askcos_site.main.views as views

# Static (not good for deployment)
# urlpatterns = static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
# urlpatterns += static(settings.MEDIA_URL,  document_root=settings.MEDIA_ROOT)

# The rest
urlpatterns = [
    # Examples:
    # re_path(r'^$', 'askcos_site.views.home', name='home'),
    # re_path(r'^blog/', include('blog.urls')),

    re_path(r'^robots.txt$', TemplateView.as_view(template_name="robots.txt", content_type="text/plain"), name="robots_file"),
    re_path(r'^examples/$', TemplateView.as_view(template_name='examples.html')),

    # Admin page
    path('admin/', admin.site.urls),

    # User pages
    path('registration/', include('django_registration.backends.one_step.urls')),
    path('registration/', include('django.contrib.auth.urls')),

    # Homepage
    re_path(r'^$', views.index, name='index'),
    re_path(r'^help/modules$', views.modules, name='modules'),
    re_path(r'^help/tutorial$', views.tutorial, name='tutorial'),
    re_path(r'^help/faq$', views.faq, name='faq'),

    # Retrosynthesis
    re_path(r'^retro/$', views.retro, name='retro_home'),
    re_path(r'^retro/target=(?P<smiles>.+)$', views.retro_target, name='retro_target'),

    # Interactive retrosynthesis
    re_path(r'^retro_interactive/$', views.retro_interactive, name='retro_interactive'),
    re_path(r'^retro_interactive/target=(?P<target>.+)$', views.retro_interactive, name='retro_interactive_target'),
    re_path(r'^retro_interactive_mcts/$', views.retro_interactive_mcts, name='retro_interactive_mcts'),
    re_path(r'^retro_interactive_mcts/target=(?P<target>.+)$', views.retro_interactive_mcts, name='retro_interactive_mcts_target'),
    re_path(r'^ajax/smiles_to_image/$', views.ajax_smiles_to_image, name='ajax_smiles_to_image'),
    re_path(r'^ajax/rxn_to_image/$', views.ajax_rxn_to_image, name='ajax_rxn_to_image'),
    re_path(r'^ajax/start_retro_mcts_celery/$', views.ajax_start_retro_mcts_celery, name='ajax_start_retro_mcts_celery'),
    re_path(r'^retro_interactive/export/(?P<_id>.+)$', views.export_retro_results, name='export_retro_results'),

    # Interactive forward prediction
    re_path(r'^synth_interactive/$', views.synth_interactive, name='synth_interactive'),

    # Site Selectivity prediction
    re_path(r'^site_prediction/$', views.site_prediction, name='site_prediction'),
    re_path(r'^ajax/get_sites/$', views.ajax_get_sites, name='ajax_get_sites'),

    # Template examination (by str(ObjectID))
    re_path(r'^template/$', views.template_view, name='template_view'),
    re_path(r'^template/target=(?P<id>.+)$', views.template_target, name='template_target'),
    re_path(r'^template/download/target=(?P<id>.+)$', views.template_target_export, name='template_target_export'),

    # Reaction examination
    # re_path(r'^reaxys/rxid=(?P<rxid>.+)$', views.rxid_target, name='rxid_target'),

    # Historians
    # re_path(ur'^history/chemicals/(?P<smiles>.+)$', views.chemical_history_check, name='chemical_history'),
    # re_path(ur'^history/reactions/(?P<smiles>.+)$', views.reaction_history_check, name='reaction_history'),

    # Pricing
    re_path(r'^buyables/$', views.buyables, name='buyables'),

    # SCScore
    re_path(r'^scscore/$', views.scscoring, name='scscoring'),
    re_path(r'^ajax/scscore_smiles/$', views.ajax_scscore_smiles, name='ajax_scscore_smiles'),

    # Drawing
    re_path(r'^draw/$', views.draw, name='draw'),
    re_path(r'^draw/smiles/(?P<smiles>.+)$', views.draw_smiles, name='draw_smiles'),
    re_path(r'^draw/template/(?P<template>.+)$', views.draw_template, name='draw_template'),
    re_path(r'^draw/reaction/(?P<smiles>.+)$', views.draw_reaction, name='draw_reaction'),
    re_path(r'^draw/mapped_reaction/(?P<smiles>.+)$', views.draw_mapped_reaction, name='draw_mapped_reaction'),
    re_path(r'^draw/highlighted_reaction/(?P<smiles>.+)$', views.draw_highlighted_reaction, name='draw_highlighted_reaction'),

    re_path(r'^draw/highlight/smiles=(?P<smiles>.+)&reacting_atoms=(?P<reacting_atoms>.+)&bonds=(?P<bonds>.+)$', views.draw_smiles_highlight, name='draw_highlight'),

    # Separation
    # re_path(r'^separation/input/$', views.sep_input, name='sep_input'),
    re_path(r'^separation/draw/(?P<fig>.+)$', views.draw_fig, name='draw_fig'),

    # Nearest Neighbor Setup
    # re_path(r'^nnRecommendation/setup/$', views.nn_predictor_setup, name='setup'),

    # Saved data
    re_path(r'^saved/$', views.user_saved_results, name='user_saved_results'),
    re_path(r'^saved/id=(?P<_id>.+)$', views.user_saved_results_id, name='user_saved_results_id'),
    re_path(r'^ajax/user_save_page/$', views.ajax_user_save_page, name='ajax_user_save_page'),
    re_path(r'^saved/delete/id=(?P<_id>.+)$', views.user_saved_results_del, name='user_saved_results_del'),

    # Blacklisted reactions
    re_path(r'^blacklisted/reactions/$', views.user_blacklisted_reactions, name='user_blacklisted_reactions'),
    re_path(r'^ajax/user_blacklist_reaction/$', views.ajax_user_blacklist_reaction, name='ajax_user_blacklist_reaction'),
    re_path(r'^blacklisted/reactions/delete/id=(?P<_id>.+)$', views.user_blacklisted_reactions_del, name='user_blacklisted_reactions_del'),
    re_path(r'^ajax/user_deactivate_reaction/$', views.ajax_user_deactivate_reaction, name='ajax_user_deactivate_reaction'),
    re_path(r'^ajax/user_activate_reaction/$', views.ajax_user_activate_reaction, name='ajax_user_activate_reaction'),

    # Blacklisted chemicals
    re_path(r'^blacklisted/chemicals/$', views.user_blacklisted_chemicals, name='user_blacklisted_chemicals'),
    re_path(r'^ajax/user_blacklist_chemical/$', views.ajax_user_blacklist_chemical, name='ajax_user_blacklist_chemical'),
    re_path(r'^blacklisted/chemicals/delete/id=(?P<_id>.+)$', views.user_blacklisted_chemicals_del, name='user_blacklisted_chemicals_del'),
    re_path(r'^ajax/user_deactivate_chemical/$', views.ajax_user_deactivate_chemical, name='ajax_user_deactivate_chemical'),
    re_path(r'^ajax/user_activate_chemical/$', views.ajax_user_activate_chemical, name='ajax_user_activate_chemical'),

    # async results
    re_path(r'^view-result/$', views.view_result, name='view_result'),
    re_path(r'^view-tree-graph/$', views.view_tree_graph, name='view_tree_graph'),
    re_path(r'^my-results/$', views.my_results, name='my_results'),

    # Reaction network
    re_path(r'^retro/network/$', views.retro_network, name='retro_network'),

    # Celery status
    re_path(r'^status/$', views.status),

    # Atom mapping
    re_path(r'^atom_mapping/$', views.atom_mapping, name='atom_mapping'),
    re_path(r'^ajax/find_atom_mapping/$', views.ajax_find_atom_mapping, name='ajax_find_atom_mapping'),

    # API endpoints
    path('api/v1/', include('askcos_site.api.urls'))
]
