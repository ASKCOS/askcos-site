from django.urls import re_path

from askcos_site import api

app_name = 'v1'

urlpatterns = [
    re_path(r'^retro/$', api.retro.singlestep, name='retro_api'),
    re_path(r'^fast-filter/$', api.fast_filter.fast_filter, name='fast_filter_api'),
    re_path(r'^context/$', api.context.neural_network, name='context_api'),
    re_path(r'^forward/$', api.forward.template_free, name='forward_api'),
    re_path(r'^impurity/$', api.impurity.impurity_predict, name='impurity_api'),
    re_path(r'^template/$', api.template.template, name='template_api'),
    re_path(r'^template/download/$', api.template.reaxys_export, name='api_template_reaxys_export'),
    re_path(r'^reactions/$', api.reactions.reactions, name='reactions_api'),
    re_path(r'^treebuilder/$', api.tree_builder.tree_builder, name='tree_builder_api'),
    re_path(r'^scscore/$', api.scscore.scscore, name='scscore_api'),
    re_path(r'^celery/$', api.status.celery_status, name='celery_api'),
    re_path(r'^celery/task/$', api.status.task_status, name='celery_task_api'),
    re_path(r'^validate-chem-name/$', api.validate_chem_name.validate_chem_name, name='validate_chem_name_api'),
    re_path(r'^buyables/search', api.buyables.buyables, name='all_buyables_api'),
    re_path(r'^buyables/add', api.buyables.add_buyable, name='add_buyables_api'),
    re_path(r'^buyables/upload', api.buyables.upload_buyable, name='upload_buyables_api'),
    re_path(r'^buyables/delete', api.buyables.delete_buyable, name='delete_buyables_api'),

    re_path(r'^cluster/$', api.cluster.cluster, name='cluster_api'),
    re_path(r'^selectivity/$', api.selectivity.selectivity, name='selectivity'),

    re_path(r'^rdkit/smiles-to-molfile/$', api.rdkit.smiles_to_molfile, name='smiles_to_molfile_api'),
    re_path(r'^rdkit/molfile-to-smiles/$', api.rdkit.molfile_to_smiles, name='molfile_to_smiles_api'),
    re_path(r'^rdkit/canonicalize/$', api.rdkit.canonicalize, name='canonicalize_api'),

    # async results
    re_path(r'^get-result/$', api.results.get_result, name='get_async_result'),
    re_path(r'^my-results/$', api.results.my_results, name='api_my_results'),
    re_path(r'^remove-result/$', api.results.remove_result, name='api_remove_results'),
    re_path(r'^poll-result/$', api.results.poll_result, name='api_poll_result'),
]