import dataclasses, json

class CustomJSONEncoder(json.JSONEncoder):
    """ Handles dataclasses. Dumps them as dicts (objects). """
    def default(self, o):
        if dataclasses.is_dataclass(o):
            return dataclasses.asdict(o)
        return super().default(o)
