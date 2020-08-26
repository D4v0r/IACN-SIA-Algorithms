from models.mini_sia import SIA, TUnit, TUnits


def main():
    sia = SIA(8, 10)
    sia.units = SIA.bone_marrow_binary(7, 8)
    print(sia.units)

    def hamming_func(s1: TUnit, s2: TUnit) -> float:
        summation = 0
        for a, b in zip(s1, s2):
            summation += a == b
        return summation

    detectors = sia.negative_selection(hamming_func, 2, 8, 2)
    print(detectors)
    antigens = ['10011111', '11101000', '10000011', '10001000', '01000001', '01011001', '11000010']
    print(sia.clonalg(hamming_func, antigens, 2, 2, 10))


if __name__ == '__main__':
    main()
