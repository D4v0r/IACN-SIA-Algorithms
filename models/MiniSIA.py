from typing import Callable, List, NewType
import numpy as np

TUnit = NewType("TUnit", str)
TUnits = List[TUnit]


class SIA:
    """Class SIA to generate TUnits and solve problems with they"""
    random_generation = lambda size: ''.join(map(str, np.random.randint(0, 2, size, int)))

    def __init__(self, g=random_generation):
        self.random_generation = g
        self.units = []

    @staticmethod
    def bone_marrow_binary(n: int, size: int) -> TUnits:
        """
        Generate n binary units of equal size.
        """
        population = ['' for _ in range(n)]
        for i in range(n):
            for num in np.random.randint(0, 2, size, int):
                population[i] += str(num)
        return population

    def negative_selection(self, affinity: Callable[[TUnit, TUnit], float], n: int, size: int,
                           threshold: int) -> TUnits:
        """
        Mature n antibodies of equal size considering the the affinity function and the threshold.
        """
        selected_units = []
        while len(selected_units) < n:
            random_unit = self.random_generation(size)
            if any(list(map(lambda x: affinity(random_unit, x) < threshold, self.units))): selected_units += [
                random_unit]
        return selected_units

    @staticmethod
    def monitoring(affinity, antibodies: TUnits, units: TUnits, threshold: int) -> TUnits:
        """
        Detect with antibodies the antigens present in a set of units considering the affinity function and the
        threshold.
        """
        antigens = set()
        for antibody in antibodies:
            for unit in units:
                if affinity(antibody, unit) >= threshold:
                    antigens.add(unit)
        return list(antigens)

    def clonalg(self, affinity: Callable[[TUnit, TUnit], float], antigens: TUnits, cloning: int, replacing: int,
                iterations: int) -> TUnits:
        """
        Memorize the antibodies useful for the antigens considering the affinity function, the percentages of cloning
        and replacing and the number of iterations.
        """
        pass
