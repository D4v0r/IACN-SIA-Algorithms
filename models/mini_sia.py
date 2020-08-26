from typing import Callable, List, NewType
import numpy as np
from random import randint

TUnit = NewType("TUnit", str)
TUnits = List[TUnit]


class SIA:
    """
    Class SIA to generate TUnits and solve problems with they

    This class is an abstraction of artificial immune systems algorithms based on the stranger model.

    ...
    Attributes
    ----------
    units : TUnits
        a collection of TUnits, represents the immune cells of this immune system.
    """

    def __init__(self, unit_size, beta):
        self.units = []
        self.memory = {}
        self.unit_size = unit_size
        self.beta = beta

    def __mutation(self, x: list, antigen, affinity) -> list:
        for i in reversed(range(len(x))):
            number_mutations = int((1 + i // affinity(x[i], antigen)) * self.beta)
            for mutation in range(number_mutations):
                rand_int = randint(0, self.unit_size - 1)
                x[i] = list(map(str, x[i]))
                x[i][rand_int] = ('1', '0')[int(x[i][rand_int])]
                x[i] = "".join(x[i])
        return x

    @staticmethod
    def __random_generation(size: int) -> TUnit:
        """"
        Generate an Random TUnit

        Args:
            size (int): size of TUnit.

        Returns:
            TUnit: new random TUnit.
        """
        return ''.join(map(str, np.random.randint(0, 2, size, int)))

    def __insert_pattern(self, antigen, new_pattern, affinity):
        best_memory = max(self.memory[antigen], new_pattern, key=lambda x: affinity(x, antigen))
        if best_memory != self.memory[antigen] and self.memory[antigen] in self.units:
            self.units.remove(self.memory[antigen])
            self.memory[antigen] = best_memory
            self.units += [best_memory]

    def __replace(self, replacing, affinity, antigen):
        random_units = self.bone_marrow_binary(replacing, self.unit_size)
        self.units = sorted(self.units, key=lambda x: affinity(x, antigen))
        for replace in range(replacing):
            self.units[replace] = random_units[replace]

    def __clonation(self, selected_antibodies):
        new_population = []
        for i in range(len(selected_antibodies)):
            number_clones = (self.beta * len(self.units)) // (i + 1)
            for iteration in range(number_clones):
                new_population += [selected_antibodies[i]]
        return new_population

    @staticmethod
    def __max_elements(antibodies, affinity, antigen, N):
        selected = []
        while len(selected) < N:
            u = max(antibodies, key=lambda x: affinity(x, antigen))
            selected += [u]
            antibodies.remove(u)
        return selected

    @staticmethod
    def bone_marrow_binary(n: int, size: int) -> TUnits:
        """
        Generate n binary units of equal size.

        Args:
            n (int): number of TUnits.
            size (int): size of each TUnit.

        Returns:
             TUnits: new TUnit population.
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

        Args:
            affinity (Callable): affinity function.
            n (int): number of detectors(antibodies) to generate.
            size (int): size of each detector(antibody).
            threshold (int): threshold to evaluate affinity of an detector.

        Returns:
             TUnits: Collection of TUnits, represents detectors(antibodies) population.
        """
        selected_units = []
        while len(selected_units) < n:
            random_unit = self.__random_generation(size)
            if any(list(map(lambda x: affinity(random_unit, x) < threshold, self.units))):
                selected_units += [random_unit]
        return selected_units

    @staticmethod
    def monitoring(affinity, antibodies: TUnits, units: TUnits, threshold: int) -> TUnits:
        """
        Detect with antibodies the antigens present in a set of units considering the affinity function and the
        threshold.

        Args:
            affinity (Callable): affinity function.
            antibodies (TUnits): represent an population of effective detectors.
            units (TUnits): represent an population of stranger units.
            threshold (int): threshold to detect an stranger.

        Returns:
            TUnits: represent a set of found pathogens.
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

        Args:
            affinity (Callable): affinity function.
            antigens (TUnits): set of dangerous units.
            cloning (int): percentage of cloning.
            replacing (int): percentage of replacing.
            iterations (int): number of cycles.

        Returns:
            TUnits
        """
        self.memory = {antigens[i]: self.units[i] for i in range(len(antigens))}

        for iteration in range(iterations):
            for antigen in antigens:
                selected_antibodies = self.__max_elements(self.units.copy(), affinity, antigen,
                                                          cloning)  # Elegir los mejores
                clones = self.__clonation(selected_antibodies)  # clonar a los mejores.
                mutated = self.__mutation(clones, antigen, affinity)  # Mutar a los que tiene mas baja afinidad.
                best_detector = max(mutated,
                                    key=lambda m: affinity(antigen, m))  # Buscar el mejor detector del antigeno.
                self.__insert_pattern(antigen, best_detector, affinity)  # Reemplazar en la memoria si es mejor.
                self.__replace(replacing, affinity, antigen)  # Reemplazar los n peores por individuos nuevos.

        return self.memory
