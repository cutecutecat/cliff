import abc


class Scenery:
    sequence: list[str]
    fitness: list[float]


class Parser(metaclass=abc.ABCMeta):
    @classmethod
    @abc.abstractmethod
    def parse(cls) -> Scenery:
        pass
