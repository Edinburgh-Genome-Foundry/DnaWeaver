import json
import networkx as nx


class JsonImportMixin:
    """Json import/export mixin for DnaSupplier."""

    @classmethod
    def from_json_data(cls, json_file=None, data=None, dna_suppliers_dict="default"):
        """

        Returns
        -------

        levels, suppliers_dict, main_id.
        """
        if dna_suppliers_dict == "default":
            dna_suppliers_dict = cls.default_suppliers_dict

        if json_file is not None:
            with open(json_file, "r") as f:
                data = json.load(f)

        if "graph" in data:
            data = data["graph"]

        sorted_suppliers, levels = _sort_suppliers(data)
        main_id = sorted_suppliers[-1]
        suppliers_dict = {}
        for supplier_id in sorted_suppliers:
            supplier_data = data[supplier_id]
            supplier_data["parameters"]["name"] = supplier_data["name"]
            supplier_data["parameters"]["suppliers"] = [
                suppliers_dict[supp_id] for supp_id in supplier_data["suppliers"]
            ]
            supplier_class = dna_suppliers_dict[supplier_data["type"]]
            supplier = supplier_class.from_dict(supplier_data["parameters"])
            supplier.id = supplier_id
            suppliers_dict[supplier_id] = supplier
        main = suppliers_dict[main_id]
        main.data = {"levels": levels, "suppliers_dict": suppliers_dict}
        return main


def _sort_suppliers(data):
    supply_graph = nx.DiGraph(
        [
            (supplier, supplier_id)
            for supplier_id, supplier_data in data.items()
            for supplier in supplier_data["suppliers"]
        ]
    )
    sorted_suppliers = []
    level = 0
    levels = {}
    while len(supply_graph):
        level += 1
        independant_suppliers = [
            n for n in supply_graph.nodes() if len(nx.ancestors(supply_graph, n)) == 0
        ]
        for supp in independant_suppliers:
            levels[supp] = level
        sorted_suppliers.extend(independant_suppliers)
        supply_graph.remove_nodes_from(independant_suppliers)
    return sorted_suppliers, levels
