"""This class implements a very thin layer to store metadata in a string."""


class SequenceString(str):
    def __new__(self, value, metadata=None):
        return str.__new__(self, value)

    def __init__(self, value, metadata=None):
        self.metadata = metadata or {}

    @staticmethod
    def from_record(record, topology="default_to_linear"):
        if topology.startswith("default_to_"):
            default_topology = topology.split("_")[-1]
            topology = record.annotations.get("topology", default_topology)
        return SequenceString(str(record.seq), metadata={"topology": topology})
