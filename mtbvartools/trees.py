import dendropy
import numpy as np
from tqdm import tqdm
from itertools import product
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


def getMeanNodeTip(ancestor_path, tree_path, remove_mask):
    tree = loadTree(tree_path)
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
    # precalculate terminal branch lengths
    print('precalculating leaf tbls....')
    with vt.CallBytestream(ancestor_path, init_calls=False, init_nodes=True) as ancestor_calls:
        leaf_tbl = {}
        for ln in tqdm(tree.leaf_nodes()):
            leaf_tbl[ln.label] = getVariantDistance(ancestor_calls, ln.parent_node.label, ln.label, remove_mask)
    print('calculating all mean tbls...')
    mean_terminal_branch_length = {}
    terminal_branch_count = {}
    for node in tqdm(list(tree.postorder_node_iter())):
        # leaf node lookup for parent node
        node_tbls = []
        for ln in node.leaf_nodes():
            node_tbls.append(leaf_tbl[ln.label])
        # calculate mean, count
        mean_terminal_branch_length[node.label] = np.mean(node_tbls)
        terminal_branch_count[node.label] = len(node_tbls)
    return mean_terminal_branch_length, terminal_branch_count


def getDescendantsPerTime(ancestor_path, tree_path, remove_mask, pseudocount=1):
    tree = loadTree(tree_path)
    # number of transmission events per time since variant node
    descendants_per_time = {}
    with vt.CallBytestream(ancestor_path, init_calls=False, init_nodes=True) as ancestor_calls:
        for node in tqdm(list(tree.postorder_internal_node_iter())):
            node_calls = ancestor_calls.nodes.loc[node.label]
            leaf_distances = []
            for leaf in node.leaf_nodes():
                leaf_calls = ancestor_calls.nodes.loc[leaf.label]
                distance = np.all([
                    node_calls != leaf_calls,
                    node_calls != 2,
                    leaf_calls != 2,
                    remove_mask != True], axis=0).sum()
                leaf_distances.append(distance)
            descendants_per_time[node.label] = len(leaf_distances) / (np.mean(leaf_distances) + pseudocount)
    return descendants_per_time


def getTransmissionClusters(ancestor_path, tree_path, remove_mask, distance):
    tree = loadTree(tree_path)
    # identify clusters where ALL strains are in a particular variant distance
    strains_in_clusters = []
    cluster_list = []
    cluster_dict = {}
    with vt.CallBytestream(ancestor_path, init_calls=False, init_nodes=True) as ancestor_calls:
        for node in tqdm(list(tree.leaf_nodes())):  # iterate across all leaf nodes
            if node.label in strains_in_clusters:
                continue  # clusters are accessible from every node, ignore already clustered strains
            cluster_set = np.asarray([node.label])
            test_node = node
            while node.parent_node is not None:
                fail_distance = False
                label_key = test_node.label  # store prior node's label to define the cluster's parent node
                test_node = test_node.parent_node  # iterate test node <- parent node
                test_set = np.asarray([n.label for n in test_node.leaf_nodes()])
                new_leaves = np.setdiff1d(test_set, cluster_set)
                comparisons = product(  # check if each new leaf is in distance of e/o and old leaves
                    new_leaves, test_set)
                for leaf1, leaf2 in comparisons:
                    if getVariantDistance(ancestor_calls, leaf1, leaf2, remove_mask) > distance:
                        fail_distance = True  # at least one strain fails distance
                        break
                if fail_distance:
                    break  # break while loop
                else:  # current set is passing
                    cluster_set = test_set
            if len(cluster_set) > 1:
                cluster_list.append(  # add cluster to flat list
                    cluster_set)
                strains_in_clusters = np.concatenate(  # add cluster strains to list of strains in SOME cluster
                    (strains_in_clusters, cluster_set))
                # store cluster keyed by defining parental node
                if label_key in cluster_dict.keys():
                    raise IndexError('Clusters should only be defined once - something is wrong!')
                cluster_dict[label_key] = cluster_set
    return cluster_dict, cluster_list, strains_in_clusters


def getChildClusterFraction(tree_path, cluster_dict):
    tree = loadTree(tree_path)
    # for each node, get the fraction of strains that are in clusters completely child to node
    clustered_fraction = {}
    defining_nodes = list(cluster_dict.keys())
    for target_node in tqdm(list(tree.postorder_internal_node_iter())):
        if target_node.label in defining_nodes:  # easy case, node defines cluster -> 100% clustered
            clustered_fraction[target_node.label] = 1
        else:  # need to do a calculation
            internal_child_nodes = [  # get all internal child nodes
                n.label for n in target_node.postorder_internal_node_iter()]
            cluster_child_nodes = np.intersect1d(  # get child nodes defining a cluster
                defining_nodes, internal_child_nodes)
            if len(cluster_child_nodes) == 0:  # no clusters found -> 0% clustered
                clustered_fraction[target_node.label] = 0
            else:
                # identify leaves that fall within a cluster
                clustered_leaves = []
                for cluster_node in cluster_child_nodes:
                    clustered_leaves += list(cluster_dict[cluster_node])
                # verify no double counts
                if len(clustered_leaves) != len(np.unique(clustered_leaves)):
                    raise ValueError('Overlapping clusters found, this will break calculation assumptions!')
                # compare to total number of child leaves
                clustered_fraction[target_node.label] = len(clustered_leaves) / len(target_node.leaf_nodes())
    return clustered_fraction