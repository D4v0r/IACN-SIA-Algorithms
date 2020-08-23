from models.mini_sia import SIA, TUnit, TUnits


def main():
    print(SIA.bone_marrow_binary(7, 4))
    sia = SIA()
    sia.units = SIA.bone_marrow_binary(7, 4)

    def hamming_func(s1: TUnit, s2: TUnit) -> float:
        summation = 0
        for a, b in zip(s1, s2):
            summation += a == b
        return summation

    detectors = sia.negative_selection(hamming_func, 2, 4, 2)
    print(detectors)


if __name__ == '__main__':
    main()
