import networkx as nx
from .dna_sources import (CommercialDnaOffer, DnaAssemblyStation, PartsLibrary,
                          PcrOutStation, DnaSourcesComparator)
DEFAULT_DNA_SOURCES_DICT = {
    'commercial':  CommercialDnaOffer,
    'assembly': DnaAssemblyStation,
    'library': PartsLibrary,
    'pcr': PcrOutStation,
    'comparator': DnaSourcesComparator,
    'main': DnaSourcesComparator
}

def _sort_suppliers(graph_data):
    supply_graph = nx.DiGraph([
        (supplier, supplier_id)
        for supplier_id, supplier_data in graph_data.items()
        for supplier in supplier_data['suppliers']
    ])
    sorted_suppliers = []
    level = 0
    levels = {}
    while len(supply_graph):
        level += 1
        independant_suppliers = [
            n for n in supply_graph.nodes()
            if len(nx.ancestors(supply_graph, n)) == 0
        ]
        for supp in independant_suppliers:
            levels[supp] = level
        sorted_suppliers.extend(independant_suppliers)
        supply_graph.remove_nodes_from(independant_suppliers)
    return sorted_suppliers, levels


def supply_network_from_json(graph_data, dna_sources_dict='default'):
    if dna_sources_dict == 'default':
        dna_sources_dict = DEFAULT_DNA_SOURCES_DICT

    sorted_suppliers, levels = _sort_suppliers(graph_data)
    main_id = sorted_suppliers[-1]
    suppliers_dict = {}
    for supplier_id in sorted_suppliers:
        supplier_data = graph_data[supplier_id]
        supplier_data['parameters']["name"] = supplier_data["name"]
        supplier_data['parameters']['suppliers'] = [
            suppliers_dict[supp_id]
            for supp_id in supplier_data['suppliers']
        ]
        supplier_class = dna_sources_dict[supplier_data["type"]]
        supplier = supplier_class.from_dict(supplier_data['parameters'])
        supplier.id = supplier_id
        suppliers_dict[supplier_id] = supplier
    return levels, suppliers_dict, main_id