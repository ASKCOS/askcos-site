# ASKCOS API v2

## Contents
- [Celery Task API](#celery-task-api)
    - [Atom mapping tool](#atom-mapping-tool)
    - [Reaction context prediction](#reaction-context-prediction)
    - [Fast filter scorer](#fast-filter-scorer)
    - [Forward prediction](#forward-prediction)
    - [Impurity prediction](#impurity-prediction)
    - [Retrosynthetic prediction](#retrosynthetic-prediction)
    - [Site selectivity prediction](#site-selectivity-prediction)
    - [General selectivity prediction](#general-selectivity-prediction)
    - [Retrosynthetic path ranking prediction](#retrosynthetic-path-ranking-prediction)
    - [Retrosynthetic tree builder tool](#retrosynthetic-tree-builder-tool)
    - [Quantum descriptor prediction](#quantum-descriptor-prediction)
- [SMILES API](#smiles-api)
    - [Canonicalize](#canonicalize)
    - [Validate](#validate)
    - [Molfile to SMILES](#molfile-to-smiles)
    - [SMILES to Molfile](#smiles-to-molfile)
- [Buyables API](#buyables-api)
    - [Main buyables endpoint](#main-buyables-endpoint)
    - [Specific buyable endpoint](#specific-buyable-endpoint)
    - [Upload buyables endpoint](#upload-buyables-endpoint)
    - [Buyables source query](#buyables-source-query)
- [Template API](#template-api)
    - [Specific template endpoint](#specific-template-endpoint)
    - [Reaxys query export](#reaxys-query-export)
    - [Template set query](#template-set-query)
- [Saved Results API](#saved-results-api)
    - [Main results endpoint](#main-results-endpoint)
    - [Specific result endpoint](#specific-result-endpoint)
    - [Check result endpoint](#check-result-endpoint)
    - [IPP result endpoint](#ipp-result-endpoint)
    - [Tree result endpoint](#tree-result-endpoint)
- [Banlist API](#banlist-api)
    - [Main banlist endpoints](#main-banlist-endpoints)
    - [Specific banlist entry endpoints](#specific-banlist-entry-endpoints)
    - [Banlist entry activation/deactivation endpoints](#banlist-entry-activationdeactivation-endpoints)
- [Authentication Token API](#authentication-token-api)
    - [Request token](#request-token)
    - [Refresh token](#refresh-token)
- [Other API](#other-api)
    - [Celery status](#celery-status)
    - [Drawing](#drawing)
    - [Reaction clustering](#reaction-clustering)
    - [Reaction lookup](#reaction-lookup)
    - [SCScorer](#scscorer)
    - [Template relevance model versions](#template-relevance-model-versions)

## Celery Task API
The API endpoints in this section are used to submit various celery tasks to the queue.
All of these endpoints require POST requests containing task parameters, and will return
the task ID associated with the created celery task.

Once the task has been created, you must retrieve the task result using the `/api/v2/celery/task` endpoint.

### Task retrieval
API endpoint for retrieving status and output of a celery task.
The task id is obtained from submission of any celery task.
Celery task results are not permanent, and will expire after 30 minutes.

URL: `/api/v2/celery/task/<task id>/`

Method: GET

Returns:

- `complete`: boolean indicating whether job is complete
- `failed`: boolean indicating if job failed
- `percent`: completion percent of job
- `message`: message regarding job status
- `state`: state of the job
- `error`: any error message if encountered
- `output`: output of celery task if complete

### Atom mapping tool
API endpoint for generating atom mappings for reactions.

URL: `/api/v2/atom-mapper/`

Method: POST

Parameters:

- `rxnsmiles` (str): reaction SMILES string
- `mapper` (str, optional): atom mapping backend to use (currently only 'WLN atom mapper')
- `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

Returns:

- `task_id`: celery task ID

### Reaction context prediction
API endpoint for context recommendation prediction using neural network model.

URL: `/api/v2/context/`

Method: POST

Parameters:

- `reactants` (str): SMILES string of reactants
- `products` (str): SMILES string of products
- `with_smiles` (bool, optional): whether to use SMILES for prediction
- `single_solvent` (bool, optional): whether to use single solvent for prediction
- `return_scores` (bool, optional): whether to also return scores
- `num_results` (int, optional): max number of results to return

Returns:

- `task_id`: celery task ID

### Fast filter scorer
API endpoint for reaction scoring using fast filter model.

URL: `/api/v2/fast-filter/`

Method: POST

Parameters:

- `reactants` (str): SMILES string of reactants
- `products` (str): SMILES string of products

Returns:

- `task_id`: celery task ID

### Forward prediction
API endpoint for template-free forward prediction task.

URL: `/api/v2/forward/`

Method: POST

Parameters:

- `reactants` (str): SMILES string of reactants
- `reagents` (str, optional): SMILES string of reagents
- `solvent` (str, optional): SMILES string of solvent
- `num_results` (int, optional): max number of results to return
- `atommap` (bool, optional): Flag to keep atom mapping from the prediction (default=False)
- `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

Returns:

- `task_id`: celery task ID

### Impurity prediction
API endpoint for impurity prediction task.

URL: `/api/v2/impurity/`

Method: POST

Parameters:

- `reactants` (str): SMILES string of reactants
- `reagents` (str, optional): SMILES string of reagents
- `products` (str, optional): SMILES string of products
- `solvent` (str, optional): SMILES string of solvent
- `top_k` (int, optional): max number of results to return
- `threshold` (float, optional): probability threshold
- `predictor` (str, optional): forward predictor to use
- `inspector` (str, optional): reaction scorer to use
- `mapper` (str, optional): reaction atom mapper to use
- `check_mapping` (bool, optional): whether to check atom mapping

Returns:

- `task_id`: celery task ID

### Retrosynthetic prediction
API endpoint for single-step retrosynthesis task.

URL: `/api/v2/retro/`

Method: POST

Parameters:

- `target` (str): SMILES string of target
- `num_templates` (int, optional): number of templates to consider
- `max_cum_prob` (float, optional): maximum cumulative probability of templates
- `filter_threshold` (float, optional): fast filter threshold
- `template_set` (str, optional): reaction template set to use
- `template_prioritizer_version` (int, optional): version number of template relevance model to use
- `precursor_prioritizer` (str, optional): name of precursor prioritizer to use (Relevanceheuristic or SCScore)
- `cluster` (bool, optional): whether or not to cluster results
- `cluster_method` (str, optional): method for clustering results
- `cluster_feature` (str, optional): which feature to use for clustering
- `cluster_fp_type` (str, optional): fingerprint type for clustering
- `cluster_fp_length` (int, optional): fingerprint length for clustering
- `cluster_fp_radius` (int, optional): fingerprint radius for clustering
- `selec_check` (bool, optional): whether or not to check for potential selectivity issues
- `attribute_filter` (list[dict], optional): template attribute filter to apply before template application
- `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

Returns:

- `task_id`: celery task ID

### Site selectivity prediction
API endpoint for site selectivity prediction task.

URL: `/api/v2/selectivity/`

Method: POST

Parameters:

- `smiles` (str): SMILES string of target

Returns:

- `task_id`: celery task ID


### General selectivity prediction
API endpoint for general selectivity prediction task.

URL: `/api/v2/gen-selectivity/`

Method: POST

Parameters:

- `reactants` (str): SMILES string of reactants
- `product` (str): SMILES string of product
- `reagents` (str, optional): SMILES string of reagents
- `solvent` (str, optional): SMILES string of solvent
- `mapped` (bool, optional): whether input is already atom mapped, default False
- `all_outcomes` (bool, optional): whether to return all outcomes, default False
- `verbose` (bool, optional): if True, return a json document, default True
- `mapper` (str, optional): which atom mapper to use ('Transformer' or 'WLN atom mapper')
- `no_map_reagents` (bool, optional): do not map reagents, default True
- `mode` (str, optional): which regioselectivity model to use ('GNN' or 'qm_GNN')

Returns:

- `task_id`: celery task ID


### Retrosynthetic path ranking prediction
API endpoint for retrosynthetic path ranking task.

URL: `/api/v2/path-ranking/`

Method: POST

Parameters:

- `trees` (str): list of trees to rank as a json string
- `cluster` (bool, optional): whether or not to cluster pathways
- `cluster_method` (str, optional): hdbscan or kmeans
- `min_samples` (int, optional): min samples for hdbscan
- `min_cluster_size` (int, optional): min cluster size for hdbscan

Returns:

- `task_id`: celery task ID


### Retrosynthetic tree builder tool
API endpoint for tree builder prediction task.

URL: `/api/v2/tree-builder/`

Method: POST

Parameters:

- `smiles` (str): SMILES string of target
- `version` (int, optional): tree builder version to use
- `max_depth` (int, optional): maximum depth of returned pathways
- `max_branching` (int, optional): maximum branching during pathway exploration
- `expansion_time` (int, optional): time limit for tree expansion
- `template_count` (int, optional): number of templates to consider
- `max_cum_prob` (float, optional): maximum cumulative probability of templates
- `buyable_logic` (str, optional): logic type for buyable termination (none/and/or)
- `max_ppg_logic` (str, optional): logic type for price based termination (none/and/or)
- `max_ppg` (int, optional): maximum price for price based termination
- `max_scscore_logic` (str, optional): logic type for synthetic complexity termination (none/and/or)
- `max_scscore` (int, optional): maximum scscore for synthetic complexity termination
- `chemical_property_logic` (str, optional): logic type for chemical property termination (none/and/or)
- `max_chemprop_c` (int, optional): maximum carbon count for termination
- `max_chemprop_n` (int, optional): maximum nitrogen count for termination
- `max_chemprop_o` (int, optional): maximum oxygen count for termination
- `max_chemprop_h` (int, optional): maximum hydrogen count for termination
- `chemical_popularity_logic` (str, optional): logic type for chemical popularity termination (none/and/or)
- `min_chempop_reactants` (int, optional): minimum reactant precedents for termination
- `min_chempop_products` (int, optional): minimum product precedents for termination
- `filter_threshold` (float, optional): fast filter threshold
- `template_set` (str, optional): template set to use
- `template_prioritizer_version` (int, optional): version number of template relevance model to use
- `buyables_source` (list[str], optional): list of source(s) to consider when looking up buyables
- `return_first` (bool, optional): whether to return upon finding the first pathway
- `max_trees` (int, optional): maximum number of pathways to return
- `score_trees` (bool, optional): whether to score trees using pathway ranking model
- `cluster_trees` (bool, optional): whether to cluster trees
- `cluster_method` (bool, optional): method to use for clustering, supports 'hdbscan' or 'kmeans'
- `cluster_min_samples` (bool, optional): minimum number of samples when using 'hdbscan'
- `cluster_min_size` (bool, optional): minimum cluster size when using 'hdbscan'
- `json_format` (str, optional): return format for trees, either 'treedata' or 'nodelink'
- `store_results` (bool, optional): whether to permanently save this result
- `description` (str, optional): description to associate with stored result
- `banned_reactions` (list[str], optional): list of reactions to not consider
- `banned_chemicals` (list[str], optional): list of molecules to not consider
- `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

Returns:

- `task_id`: celery task ID

### Quantum descriptor prediction
API endpoint for descriptor-predictor prediction task.

URL: `/api/v2/descriptors/`

Method: POST

Parameters:

- `smiles` (str): SMILES string

Returns:

- `task_id`: celery task ID


## SMILES API
The API endpoints in this section provide various utilities for working with SMILES.

### Canonicalize
Canonicalize the specified SMILES using RDKit.

URL: `/api/v2/rdkit/smiles/canonicalize/`

Method: POST

Parameters:

- `smiles` (str): SMILES string to canonicalize
- `isomericSmiles` (bool, optional): whether to generate isomeric SMILES

Returns:

- `smiles` (str): canonicalized SMILES

### Validate
Check the syntax and validity of a SMILES string.

URL: `/api/v2/rdkit/smiles/validate/`

Method: POST

Parameters:

- `smiles` (str): SMILES string

Returns:

- `correct_syntax` (bool): correctness of the SMILES syntax
- `valid_chem_name` (bool): validity of the chemical

RDKit ref: https://github.com/rdkit/rdkit/issues/2430

### Molfile to SMILES
Convert the provided Molfile to a SMILES string.

URL: `/api/v2/rdkit/smiles/from_molfile/`

Method: POST

Parameters:

- `molfile` (str): Molfile input
- `isomericSmiles` (bool, optional): whether to generate isomeric SMILES

Returns:

- `smiles` (str): canonical SMILES

### SMILES to Molfile
Convert the provided SMILES string to a Molfile.

URL: `/api/v2/rdkit/smiles/to_molfile/`

Method: POST

Parameters:

- `smiles` (str): SMILES input

Returns:

- `molfile` (str): Molfile output


## Buyables API
The API endpoints in this section provide methods for accessing and modifying the buyables database.
User authentication is required for any endpoints which modify the database.

### Main buyables endpoint
API endpoint for viewing buyables or creating new buyables.

URL: `/api/v2/buyables/`

Method: GET

Query Parameters:

- `q` (str, optional): search query, e.g. SMILES string
- `source` (list, optional): list of source(s) to consider when looking up buyables
- `regex` (bool, optional): whether or not to treat `q` as regex pattern (default: False)
- `returnLimit` (int, optional): maximum number of results to return (default: 100)
- `canonicalize` (bool, optional): whether or not to canonicalize `q` (default: True)

Returns:

- `search`: query pattern used for search
- `result`: list of buyables matching search query

Method: POST

Parameters:

- `smiles` (str): SMILES string of buyable
- `ppg` (float): price of buyable
- `source` (str, optional): source of data (default: '')
- `allowOverwrite` (bool, optional): whether or not to overwrite existing duplicates (default: True)

Returns:

- `success`: true if buyable was created successfully
- `error`: error message if not successful
- `duplicate`: true if buyable was not added because it already existed
- `inserted`: buyables entry which was created if it didn't exist
- `updated`: buyables entry which was updated if it existed and `allowOverwrite = True`

### Specific buyable endpoint
API endpoint for viewing or deleting specific buyable.
The buyable ID can be obtained from querying the main buyables endpoint.

URL: `/api/v2/buyables/<buyable id>/`

Method: GET

Returns:

- `_id`: the requested buyable id
- `result`: the requested buyable
- `error`: error message if encountered

Method: DELETE

Returns:

- `success`: true if deletion was successful
- `error`: error message if encountered

### Upload buyables endpoint
API endpoint for uploading buyables data.

URL: `/api/v2/buyables/upload/`

Method: POST

Parameters:

- `upload_file` (file): file containing buyables data
- `format` (str): file format, either json or csv
- `returnLimit` (int): maximum number of results to return
- `allowOverwrite` (bool): whether or not to overwrite existing duplicates

Returns:

- `success`: true if buyable was created successfully
- `error`: error message if not successful
- `inserted`: list of new buyable entries, up to `returnLimit`
- `updated`: list of updated buyable entries, up to `returnLimit`
- `inserted_count`: total number of inserted entries
- `updated_count`: total number of updated entries
- `duplicate_count`: total number of duplicate entries if `allowOverwrite = False`
- `count`: total number of successfully uploaded entries
- `total`: total number of uploaded entries, including errors

### Buyables source query
API endpoint to query available buyables sources.

URL: `/api/v2/buyables/sources/`

Method: GET

Returns:

- `sources`: list of available buyables sources


## Template API
The API endpoints in this section are for accessing retrosynthetic template data.

### Specific template endpoint
API endpoint for accessing a particular template entry.
The template ID is returned from certain tasks, like the one-step retrosynthetic prediction.

URL: `/api/v2/template/<template id>/`

Method: GET

Returns:

- `template`: reaction template

### Reaxys query export
API endpoint for exporting template reference data as a Reaxys query.
Note that this only applies for the `reaxys` template set.

URL: `/api/v2/template/<template id>/export/`

Method: GET

Returns:

- JSON format Reaxys query

### Template set query
API endpoint to query available template sets.

URL: `/api/v2/template/sets/`

Method: GET

Returns:

- `template_sets`: list of available template sets
- `attributes`: list of attribute names for each available template set

## Saved Results API
The API endpoints in this section are for accessing saved results.
Note that these are different from the celery task results.
Currently, saved results include async tree builder jobs submitted via the web client
and any web pages which have been manually saved.
User authentication is required to access these endpoints.

### Main results endpoint
API endpoint for accessing a user's job results.

URL: `/api/v2/results/`

Method: GET

Returns:

- `results`: list of results belonging to the current user

### Specific result endpoint
API endpoint for accessing or deleting a particular result.
The result ID can be obtained from the main results endpoint.
 
URL: `/api/v2/results/<result id>/`

Method: GET

Returns:

- `id`: the requested result id
- `result`: the requested result
- `error`: error message if encountered

Method: DELETE

Returns:

- `success`: true if deletion was successful
- `error`: error message if encountered

### Check result endpoint
API endpoint for checking the status of a particular result.
The result ID can be obtained from the main results endpoint.

URL: `/api/v2/results/<result id>/check/`

Method: GET

Returns:

- `state`: current state of the job
- `error`: error message if encountered

### IPP result endpoint
API endpoint which processes results for display in IPP.
Returns a single tree obtained by merging individual pathways and a dictionary of precursor results.
The result ID can be obtained from the main results endpoint.

URL: `/api/v2/results/<result id>/ipp/`

Method: GET

Returns:

- `id`: the requested result id
- `result`: the requested result, containing `settings` and `result`
- `error`: error message if encountered

### Tree result endpoint
API endpoint which processes results for display in tree results page.
Returns list of trees in nodelink format.
The result ID can be obtained from the main results endpoint.

URL: `/api/v2/results/<result id>/tree/`

Method: GET

Returns:

- `id`: the requested result id
- `result`: the requested result, containing `settings` and `result`
- `error`: error message if encountered


## Banlist API
The API endpoints in this section are for accessing and modifying user chemical and reaction banlists.
User authentication is required to access these endpoints.

### Main banlist endpoints

URLs:

- `/api/v2/banlist/chemicals`
- `/api/v2/banlist/reactions`

Method: GET

Returns: list of banned chemical or reaction entries belonging to the currently authenticated user

Method: POST

Parameters:

- `smiles` (str): chemical or reaction SMILES string to be added
- `description` (str, optional): text description
- `datetime` (str, optional): timestamp, default current time
- `active` (bool, optional): whether this entry is active, default true

Returns: created entry

### Specific banlist entry endpoints
API endpoints for accessing or deleting a particular banlist entry.
The entry ID can be obtained from the main banlist endpoints.

URLs:

- `/api/v2/banlist/chemicals/<id>`
- `/api/v2/banlist/reactions/<id>`

Method: GET

Returns: entry with requested id

Method: DELETE

Returns:

- `success`: true if successfully deleted
- `data`: data from deleted entry

### Banlist entry activation/deactivation endpoints
API endpoints for activating or deactivating a particular banlist entry.
The entry ID can be obtained from the main banlist endpoints.

URLs:

- `/api/v2/banlist/chemicals/<id>/activate`
- `/api/v2/banlist/chemicals/<id>/deactivate`
- `/api/v2/banlist/reactions/<id>/activate`
- `/api/v2/banlist/reactions/<id>/deactivate`


Method: GET

Returns:

- `success`: true if successfully activated
- `data`: updated entry data


## Authentication Token API
The API endpoints in this section are for obtaining JSON Web Tokens to use for authenticating
to API endpoints which require authentication.

### Request token
API endpoint for requesting a temporary authentication token.

URL: `/api/v2/token-auth/`

Method: POST

Parameters:

- `username`: username
- `password`: password

Returns:
- `token`: JSON Web Token that can be used for authenticated requests

### Refresh token
API endpoint for obtaining a refreshed token. Useful for extending expiration time.

URL: `/api/v2/token-refresh/`

Method: POST

Parameters:

- `token`: valid, non-expired token

Returns:
- `token`: JSON Web Token that can be used for authenticated requests


## Other API

### Celery status
API endpoint for retrieving celery worker status.

URL: `/api/v2/celery/`

Method: GET

Returns:

- `queues`: list of worker information for each celery queue

### Drawing
API endpoint for drawing molecules, reactions, and templates.

Notes:

- Both GET and POST requests are possible. GET requests may be easier
  for simple linking, while POST requests are better for complex data.
- If `input_type` is not specified, will attempt to determine type.
  Specifying `input_type` can provide faster results.

URL: `/api/v2/draw/`

Method: GET, POST

Parameters:

- `smiles` (str): input SMILES (or SMARTS) string
- `input_type` (str, optional): one of 'chemical', 'reaction', or 'template'
- `transparent` (bool, optional): whether background should be transparent (chemical only)
- `draw_map` (bool, optional): whether atom mapping should be drawn (reaction only)
- `highlight` (bool, optional): whether to highlight mapped atoms (reaction or chemical)
- `reacting_atoms` (list, optional): list of atom scores to highlight and label (chemical only)

Returns: PNG image of input SMILES

### Reaction clustering
API endpoint for clustering similar transformed outcomes

URL: `/api/v2/cluster/`

Method: POST

Parameters:

- `original` (str): smiles string
- `outcomes` (list): list of smiles strings of outcomes
- `feature` (str, optional): features to use [original', 'outcomes', 'all']
- `fingerprint` (str, optional): fingerprint type ['morgan']
- `fpradius` (int, optional): fingerprint radius, default 1
- `fpnbits` (int, optional): fingerprint bits, default 512
- `clustermethod` (str, optional): cluster method ['hdbscan', 'kmeans']
- `scores` (list, optional): list of scores of precursors

Returns:

- `request`: dictionary of request parameters
- `output`: list of cluster indices for outcomes

### Reaction lookup
API endpoint for reaction lookup task.

URL: `/api/v2/reactions/`

Method: POST

Parameters:

- `ids` (list): list of reaction ids to retrieve
- `template_set` (str, optional): template set to search within

Returns:

- `reactions`: list of reactions

### SCScorer
API endpoint for scscore prediction task.

URL: `/api/v2/scscore/`

Method: POST

Parameters:

- `smiles` (str): SMILES string of target

Returns:

- `output`: synthetic complexity score

### Template relevance model versions
API endpoint for querying available retrosynthetic models for a given template set.

URL: `/api/v2/retro/models/`

Method: GET

Parameters:

- `template_set` (str): template set name

Returns:

- `versions`: List of version numbers that are available
