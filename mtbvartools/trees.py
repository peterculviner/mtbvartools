import dendropy

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
    taxon_labels = []
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
