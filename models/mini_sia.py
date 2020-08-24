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

    def __init__(self):
        self.units = []
        self.memory = []

        
    @staticmethod
    def __mutation(x: list,size: int) -> list:
        #x = list(map(str,x))
        for i in range(len(x)):
            rand_int=randint(0,size-1)
            x[i] = list(map(str,x[i]))
            x[i][rand_int] = ('1','0')[int(x[i][rand_int])]
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

    @staticmethod
    def __max_elements(list0, list1, N):
        final_list = [] 
    
        for i in range(0, N):
            max1 = 0
            for j in range(len(list1)):
                if list1[j] > max1:
                    max1 = list1[j]                    
                    ele1 = list0[j]
            list1.remove(max1)            
            list0.remove(ele1)
            final_list.append(ele1)
        return final_list

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
        memory = {x:self.units[0] for x in antigens}
        
        for iteration in range(iterations):
            for antigen in antigens:
                aff = list(map(lambda x: affinity(antigen, x) , self.units))                
                C = self.__max_elements(self.units.copy(), aff, cloning)
                muted = self.__mutation(C,len(self.units[0]))                
                aff = list(map(lambda x: affinity(antigen, x) , muted))
                max_affinity = max(muted,key=lambda m: affinity(antigen, m))
                memory[antigen] = max(memory[antigen],max_affinity)
        return memory
                
                
                
                
                
                
        
                
                
                


        return
