import dendropy
import numpy as np
from tqdm import tqdm
import mtbvartools as vt

def writeLabelledTree(input_path=None, output_path=None, input_tree=None, **kwargs):
    """
    Label internal nodes of input_path (or input_tree) in postorder as internal<n>.
    Returns tree. If output path provided, also writes tree to newick - dendropy write func args:
        suppress_leaf_taxon_labels=False,
        suppress_leaf_node_labels=True,
        suppress_internal_taxon_labels=False,
        suppress_internal_node_labels=False,
        suppress_rooting=False,
        suppress_edge_lengths=False,
        unquoted_underscores=False,
        preserve_spaces=False,
        store_tree_weights=False,
        taxon_token_map=None,
        suppress_annotations=True,
        annotations_as_nhx=False,
        suppress_item_comments=True,
        node_label_element_separator=' ',
        node_label_compose_fn=None,
        edge_label_compose_fn=None,
        real_value_format_specifier='',
        ignore_unrecognized_keyword_arguments=False
    """
    if input_tree is None and input_path is not None:
        tree = dendropy.Tree.get(
            path=input_path, schema='newick', preserve_underscores=True)
    elif input_tree is not None and input_path is None:
        tree = input_tree
    else:
        raise ValueError('Exactly one of "input_path" and "input_tree" is required as input.')
    i = 1
    for node in tree.postorder_node_iter():
        if node.taxon == None:
            node.taxon = dendropy.Taxon(label=f'internal{i}')
            tree.taxon_namespace.add_taxon(node.taxon)
            i += 1
    if output_path is not None:
        tree.write(
            path=output_path,
            schema='newick',
            **kwargs)
    return tree


def getVariantDistance(ancestor_calls, node1, node2, remove_mask):
    n1_calls = ancestor_calls.nodes.loc[node1]
    n2_calls = ancestor_calls.nodes.loc[node2]
    distance = np.all([
        n1_calls != n2_calls,
        n1_calls != 2,
        n2_calls != 2,
        remove_mask != True], axis=0).sum()
    return distance


def loadTree(nwk_path):
    # load tree
    tree = dendropy.Tree.get(
        path=nwk_path,
        schema='newick', preserve_underscores=True)
    # copy leaf taxon labels -> leaf label
    for leaf in tree.leaf_nodes():
        leaf.label = leaf.taxon.label
    return tree


def getMeanNodeTip(ancestor_path, tree_path, remove_mask, allow_terminal=False):
    tree = loadTree(tree_path)
    if allow_terminal:
        # preload mean node tip allowing terminal branches
        mean_node_tip = {}
        with vt.CallBytestream(ancestor_path, init_calls=False, init_nodes=True) as ancestor_calls:
            for node in tqdm(list(tree.postorder_node_iter())):
                if node.is_internal():
                    if node.parent_node == None:
                        node_tip_distances = []  # root node has no parent
                    else:
                        node_tip_distances = [np.mean([  # calculate initial value of parental node to child leaves
                                getVariantDistance(ancestor_calls, node.parent_node.label, leaf_node.label, remove_mask)
                                for leaf_node in node.leaf_nodes()])]
                    for n in node.postorder_iter():
                        if n != node:
                            node_tip_distances.append(  # add precalculated child node values to the list
                                mean_node_tip[n.label])
                    mean_node_tip[node.label] = np.mean(  # finally calculate the value at this node
                        node_tip_distances)
                else:
                    mean_node_tip[node.label] = getVariantDistance(ancestor_calls, node.parent_node.label, node.label, remove_mask)
        return mean_node_tip
    else:
        # preload mean node tip, excluding terminal branches
        mean_node_tip_internal = {}
        with vt.CallBytestream(ancestor_path, init_calls=False, init_nodes=True) as ancestor_calls:
            for node in tqdm(list(tree.postorder_internal_node_iter())):
                node_tip_distances = [np.mean([  # calculate initial value of parental node to child leaves
                    getVariantDistance(ancestor_calls, node.label, leaf_node.label, remove_mask)
                    for leaf_node in node.leaf_nodes()])]
                for n in node.postorder_internal_node_iter():
                    if n != node:
                        node_tip_distances.append(  # add precalculated child node values to the list
                            mean_node_tip_internal[n.label])
                mean_node_tip_internal[node.label] = np.mean(  # finally calculate the value at this node
                    node_tip_distances)
        return mean_node_tip_internal


def getMeanTerminalBranchLengths(ancestor_path, tree_path, remove_mask):
    tree = loadTree(tree_path)
    # mean terminal branch length
    mean_terminal_branch_length = {}
    with vt.CallBytestream(ancestor_path, init_calls=False, init_nodes=True) as ancestor_calls:
        for node in tqdm(list(tree.postorder_node_iter())):
            leaf_node_tbls = []
            for ln in node.leaf_nodes():
                leaf_node_tbls.append(
                    getVariantDistance(ancestor_calls, ln.parent_node.label, ln.label, remove_mask))
            mean_terminal_branch_length[node.label] = np.mean(
                leaf_node_tbls)
    return mean_terminal_branch_length