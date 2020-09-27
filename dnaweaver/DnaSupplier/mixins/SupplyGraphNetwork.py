class SupplyNetworkMixin:
    def compute_supply_graph(self):
        """Return elements to plot the supply graph underlying this DnaSupplier.

        Returns
        -------

        edges
          A list [(s1,s2), (s1,s3), (s2, s5)...] of couples of DnaSuppliers in a
          supplier-supplied relationship.

        levels
          A list of lists [[s1,s2], [s4,s8,s9]...] of sources. The first
          sublist (first level) are all sources at the farthest distance from
          the current source in the supply graph, and the last sublist contains
          only the current DnaSupplier.
        """

        source_max_level = {}
        edges = []

        def rec(source, depth, seen_sources):
            if source in seen_sources:
                return
            if source not in source_max_level:
                source_max_level[source] = depth
            else:
                source_max_level[source] = max(source_max_level[source], depth)
            new_seen_sources = seen_sources + [source]
            if hasattr(source, "suppliers"):
                for other in source.suppliers:
                    edges.append((other, source))
                    rec(other, depth + 1, new_seen_sources)
            elif hasattr(source, "supplier"):
                edges.append((source.supplier, source))
                rec(source.supplier, depth + 1, new_seen_sources)
            if hasattr(source, "primers_supplier"):
                edges.append((source.primers_supplier, source))
                rec(source.primers_supplier, depth + 1, new_seen_sources)

        rec(self, depth=0, seen_sources=[])
        levels = [
            [source for source, level in source_max_level.items() if level == i]
            for i in range(max(source_max_level.values()) + 1)
        ][::-1]

        return edges, levels

    def dict_supply_graph(self):
        sources = {}

        def rec(source, depth=0):

            if source in sources:
                return
            if hasattr(source, "is_ghost_source") and source != self:
                return
            sources[source.name] = source.dict_description()
            sources[source.name]["_depth"] = depth

            providers = sources[source.name]["providers"] = []
            if hasattr(source, "suppliers"):
                for other in source.suppliers:
                    providers.append(other.name)
                    rec(other, depth + 1)
            # if hasattr(source, "dna_supplier"):
            #     providers.append(source.dna_supplier.name)
            #     rec(source.dna_supplier, depth + 1)
            # if hasattr(source, "primers_supplier"):
            #     providers.append(source.primers_supplier.name)
            #     rec(source.primers_supplier, depth + 1)
            # if hasattr(source, "dna_suppliers"):
            #     for other in source.dna_suppliers:
            #         providers.append(other.name)
            #         rec(other, depth + 1)

        rec(self)
        return sources
