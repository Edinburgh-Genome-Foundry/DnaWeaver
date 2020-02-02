from copy import deepcopy
from .ObjectDict import ObjectDict
from . import mixins


class AssemblyPlanReport(
    mixins.PlotsMixin,
    mixins.FolderReportMixin,
    mixins.GenbankExportMixin,
    mixins.PdfReportMixin,
):
    def __init__(self, plan, sources):
        self.plan = ObjectDict.from_dict(plan)
        self.sources = ObjectDict.from_dict(sources)

    @staticmethod
    def from_dnaweaver_quote(quote):
        plan = quote.assembly_plan_as_dict()
        sources = quote.source.dict_supply_graph()
        return AssemblyPlan(plan, sources)

    def to_steps_list(self):
        plan = deepcopy(self.plan)
        nodes = []

        def rec(node, depth=0):
            if node.get("_visited", False):
                return
            node["_visited"] = True
            assembly_plan = node.get("assembly_plan", [])
            node["children"] = [n["id"] for n in assembly_plan]
            nodes.append(node)
            for other in sorted(
                assembly_plan, key=lambda n: n["segment_start"]
            ):
                rec(other)

        rec(plan)
        return nodes
