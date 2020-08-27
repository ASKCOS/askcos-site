import networkx as nx


def graph_to_results(data):
    """
    Convert new tree builder data structure to IPP results.

    Args:
        data (dict): graph in networkx node link json format

    Returns:
        dict: keys: chemical smiles
              values: list of dictionaries containing reaction details
    """
    graph = nx.node_link_graph(data)

    results = {}
    for node, node_data in graph.nodes.items():
        if node_data['type'] == 'reaction':
            continue

        precursors = []
        for rxn in graph.successors(node):
            rxn_data = graph.nodes[rxn]
            precursors.append({
                'rank': rxn_data['rank'],
                'smiles_split': rxn_data['precursor_smiles'].split('.'),
                'smiles': rxn_data['precursor_smiles'],
                'plausibility': rxn_data['plausibility'],
                'template_score': rxn_data['template_score'],
                'score': rxn_data['template_score'],
                'templates': rxn_data['tforms'],
                'num_examples': rxn_data['num_examples'],
                'necessary_reagent': rxn_data['necessary_reagent'],
                'rms_molwt': rxn_data['rms_molwt'],
                'num_rings': rxn_data['num_rings'],
                'scscore': rxn_data['scscore'],
            })

        if not precursors:
            continue

        precursors.sort(key=lambda x: x['rank'])

        results[node] = precursors

    return results


def combine_trees(trees):
    """
    Combine a list of retrosynthetic trees into a unified tree.

    Args:
        trees (list): list of pathways in networkx tree data json format

    Returns:
        graph in networkx node link json format
    """
    graph = nx.compose_all([nx.tree_graph(tree) for tree in trees])

    for node, node_data in graph.nodes.items():
        if node_data.get('is_chemical'):
            node_data['type'] = 'chemical'
            del node_data['is_chemical']
        elif node_data.get('is_reaction'):
            node_data['type'] = 'reaction'
            del node_data['is_reaction']

    attrs = {'source': 'from', 'target': 'to', 'name': 'id', 'key': 'key', 'link': 'edges'}
    return nx.node_link_data(graph, attrs=attrs)
