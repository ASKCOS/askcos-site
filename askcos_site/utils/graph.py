import networkx as nx

from askcos.retrosynthetic.mcts.utils import generate_unique_node
from askcos.utilities.descriptors import rms_molecular_weight, number_of_rings
from askcos_site.globals import scscorer, retro_transformer

NODE_LINK_ATTRS = {'source': 'from', 'target': 'to', 'name': 'id', 'key': 'key', 'link': 'edges'}


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
    if not trees:
        return

    graph = nx.compose_all([nx.tree_graph(tree) for tree in trees])

    for node, node_data in graph.nodes.items():
        if node_data.get('is_chemical'):
            node_data['type'] = 'chemical'
            del node_data['is_chemical']
        elif node_data.get('is_reaction'):
            node_data['type'] = 'reaction'
            del node_data['is_reaction']

    return nx.node_link_data(graph, attrs=NODE_LINK_ATTRS)


def chem_to_results(data, tree=None, template_set=None):
    """
    Convert old tree builder data structure to IPP results.

    If a tree is provided (in networkx node link format), reaction nodes in the
    tree are updated with the precursor rank based on the results in ``data``.

    Args:
        data (dict): keys: chemical smiles
                     values: list of dictionaries containing reaction details
        tree (dict, optional): graph in networkx node link format
        template_set (str, optional): template set for retrieving metadata

    Returns:
        dict: keys: chemical smiles
              values: list of dictionaries containing reaction details
    """
    results = {}
    for smiles, reactions in data.items():
        # Sort results from high score to low score
        reactions.sort(key=lambda x: x['template_score'], reverse=True)

        precursors = {}
        for i, rxn in enumerate(reactions):
            precursor_smiles = '.'.join(rxn['reactant_smiles'])
            if precursor_smiles in precursors:
                continue

            template_data = None
            if tree is not None:
                rxn_smiles = precursor_smiles + '>>' + smiles
                for node in tree['nodes']:
                    if node['smiles'] == rxn_smiles:
                        # Add rank attribute to the tree
                        node['rank'] = i + 1
                        # Retrieve template data
                        template_data = {
                            'tforms': node['tforms'],
                            'num_examples': node['num_examples'],
                            'necessary_reagent': node['necessary_reagent'],
                        }

            if template_data is None:
                template_data = retro_transformer.retrieve_template_metadata(rxn['tforms'], template_set=template_set)

            precursors[precursor_smiles] = {
                'rank': i + 1,
                'smiles_split': rxn['reactant_smiles'],
                'smiles': precursor_smiles,
                'plausibility': rxn['plausibility'],
                'template_score': rxn['template_score'],
                'score': rxn['template_score'],
                'templates': template_data['tforms'],
                'num_examples': template_data['num_examples'],
                'necessary_reagent': template_data['necessary_reagent'],
                'rms_molwt': rms_molecular_weight(precursor_smiles),
                'num_rings': number_of_rings(precursor_smiles),
                'scscore': scscorer.get_score_from_smiles(precursor_smiles, noprice=True),
            }

        results[smiles] = list(precursors.values())

    return results


def combine_old_trees(trees):
    """
    Combine a list of retrosynthetic trees into a unified tree. Node IDs are
    reassigned to ensure the result is a prefix tree.

    Args:
        trees (list): list of pathways in networkx tree data json format

    Returns:
        graph in networkx node link json format
    """
    if not trees:
        return

    def _helper(_trees, _root, _graph):
        if not _trees:
            return
        nodes = {}
        children = {}
        for tree in _trees:
            # Make a shallow copy of the tree
            node = dict(tree)
            # Remove the children attribute
            rest = node.pop('children')
            # This stores the node attributes
            nodes[node['smiles']] = node
            # This stores the children of this node
            children.setdefault(node['smiles'], []).extend(rest)
        for smiles, attributes in nodes.items():
            new_id = generate_unique_node()
            _graph.add_node(new_id, **attributes)
            _graph.add_edge(_root, new_id)
            _helper(children[smiles], new_id, _graph)

    graph = nx.DiGraph()

    # Add a temporary root node to help construct the tree
    root = nx.utils.generate_unique_node()
    graph.add_node(root)

    # Merge the individual trees into a prefix tree
    _helper(trees, root, graph)

    # Relabel root with NIL UUID so we can easily identify it
    real_root = next(graph.successors(root))
    nx.relabel_nodes(graph, {real_root: '00000000-0000-0000-0000-000000000000'}, copy=False)

    graph.remove_node(root)

    for node, node_data in graph.nodes.items():
        if node_data.get('is_chemical'):
            node_data['type'] = 'chemical'
            del node_data['is_chemical']
        elif node_data.get('is_reaction'):
            node_data['type'] = 'reaction'
            del node_data['is_reaction']

    return nx.node_link_data(graph, attrs=NODE_LINK_ATTRS)
