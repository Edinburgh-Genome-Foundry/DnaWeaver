from collections import defaultdict
from .DnaSource import (ExternalDnaOffer, DnaAssemblyStation,
                        PcrOutStation, PartsLibrary)
from .plotting import plot_assembly_tree
import pandas as pd


def tree_to_assembly_protocol(tree):

    sources_counter = defaultdict(lambda: 0)
    counter = [1]
    steps = []

    def minify_name(source):
        return source.name.replace(" ", "_").replace(".", "_").lower()

    def rec(node):

        quote, children = node
        sources_counter[quote.source] += 1
        construct_name = (minify_name(quote.source) +
                          "%04d" % sources_counter[quote.source])
        children_fragments = [
            rec(child)
            for segment, child in children.items()
        ]

        steps.append({
            "TaskNumber": counter[0],
            "SequenceName": construct_name,
            "Type": {ExternalDnaOffer: "Order",
                     DnaAssemblyStation: "Assembly",
                     PcrOutStation: "PCR-out",
                     PartsLibrary: "Library"}[quote.source.__class__],
            "Source": quote.source.name,
            "Sequence/fragments": (quote.sequence if children_fragments == []
                                   else ",".join(children_fragments))
        })
        counter[0] += 1
        return construct_name
    rec(tree)

    return pd.DataFrame.from_records(
        steps,
        columns=["TaskNumber", "SequenceName", "Type",
                 "Source", "Sequence/fragments"]
    )

def assembly_protocol_to_spreadsheet(assembly_protocol, filename):
    orders = assembly_protocol[assembly_protocol.Type == "Order"].sort(
        columns=["Source", "TaskNumber"]
    )
    assemblies = assembly_protocol[assembly_protocol.Type == "Assembly"].sort(
        columns=["TaskNumber"]
    )
    with open(filename, "w+") as f:
        f.write("ASSEMBLY PLANNING\n\nAssemblies\n\n")
        assemblies.to_csv(
            f,
            columns=["Source", "SequenceName", "Sequence/fragments"],
            index=False
        )
        f.write("\n\nOrders\n\n")
        orders.to_csv(
            f,
            columns=["Source", "SequenceName", "Sequence/fragments"],
            index=False
        )
