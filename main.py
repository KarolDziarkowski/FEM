import numpy as np
import math as mat


class Elem4:
    l_elem = 4
    x = np.array([0.0, 0.025, 0.025, 0.0])
    y = np.array([0.0, 0.0, 0.025, 0.025])

    n = 1 / mat.sqrt(3)
    pcE = np.array([-n, n, n, -n])
    pcN = np.array([-n, -n, n, n])

    dNdN = np.zeros((4, 4))
    dNdE = np.zeros((4, 4))
    dNdX = np.zeros((4, 4))
    dNdY = np.zeros((4, 4))

    detJ = np.zeros(4)

    jakob = np.zeros((4, 2, 2))
    jakob_odwr = np.zeros((4, 2, 2))
    suma_H = np.zeros((4, 4))

    @staticmethod
    def pochodne():
        for x in range(Elem4.l_elem):
            wsp = 0.25
            Elem4.dNdE[x, 0] = wsp * (1 - Elem4.pcN[x])
            Elem4.dNdE[x, 1] = -wsp * (1 - Elem4.pcN[x])
            Elem4.dNdE[x, 2] = -wsp * (1 + Elem4.pcN[x])
            Elem4.dNdE[x, 3] = wsp * (1 + Elem4.pcN[x])

            Elem4.dNdN[x, 0] = wsp * (1 - Elem4.pcE[x])
            Elem4.dNdN[x, 1] = wsp * (1 + Elem4.pcE[x])
            Elem4.dNdN[x, 2] = -wsp * (1 + Elem4.pcE[x])
            Elem4.dNdN[x, 3] = -wsp * (1 - Elem4.pcE[x])

        print("Pochodne funkcji kształtu dNdE:\n", Elem4.dNdE)
        print("Pochodne funkcji kształtu dNdN:\n", Elem4.dNdN)

    @staticmethod
    def jakobian():
        dXdE = 0.0
        dYdE = 0.0
        dXdN = 0.0
        dYdN = 0.0

        for x in range(Elem4.l_elem):
            for y in range(4):
                dXdE += Elem4.dNdE[x, y] * -Elem4.x[y]
                dXdN += Elem4.dNdN[x, y] * -Elem4.x[y]
                dYdE += Elem4.dNdE[x, y] * -Elem4.y[y]
                dYdN += Elem4.dNdN[x, y] * -Elem4.y[y]
            Elem4.jakob[x, 0, 0] = dXdE
            Elem4.jakob[x, 0, 1] = dXdN
            Elem4.jakob[x, 1, 0] = dYdE
            Elem4.jakob[x, 1, 1] = dYdN
            dXdE = 0.0
            dXdN = 0.0
            dYdE = 0.0
            dYdN = 0.0

        for x in range(Elem4.l_elem):
            Elem4.jakob_odwr[x, 0, 0] = Elem4.jakob[x, 1, 1]
            Elem4.jakob_odwr[x, 0, 1] = Elem4.jakob[x, 1, 0]
            Elem4.jakob_odwr[x, 1, 0] = Elem4.jakob[x, 0, 1]
            Elem4.jakob_odwr[x, 1, 1] = Elem4.jakob[x, 0, 0]

        for x in range(Elem4.l_elem):
            Elem4.detJ[x] = Elem4.jakob[x, 0, 0] * Elem4.jakob[x, 1, 1] - Elem4.jakob[x, 0, 1] * Elem4.jakob[x, 1, 0]

        print("\nJakobian:\n", Elem4.jakob)


    @staticmethod
    def pochodne2():
        for x in range(Elem4.l_elem):
            for y in range(4):
                Elem4.dNdX[x, y] = -(1 / Elem4.detJ[x]) * (
                        Elem4.jakob_odwr[x, 0, 0] * Elem4.dNdE[x, y] + Elem4.jakob_odwr[x, 0, 1] * Elem4.dNdN[x, y])
                Elem4.dNdY[x, y] = -(1 / Elem4.detJ[x]) * (
                        Elem4.jakob_odwr[x, 1, 0] * Elem4.dNdE[x, y] + Elem4.jakob_odwr[x, 1, 1] * Elem4.dNdN[x, y])
        print("dndx\n", Elem4.dNdX)
        print("dndy\n", Elem4.dNdY)

    @staticmethod
    def macierzH():
        H = np.zeros((4, 4, 4))
        wsp_K = 30.0
        for x in range(Elem4.l_elem):
            for y in range(4):
                for z in range(4):
                    H[x, y, z] = wsp_K * (Elem4.dNdX[x, y] * Elem4.dNdX[x, z] + Elem4.dNdY[x, y] * Elem4.dNdY[x, z]) * Elem4.detJ[x]

        print("\nMacierze H przed sumowaniem:\n", H)

        for x in range(Elem4.l_elem):
            for y in range(4):
                Elem4.suma_H[x, y] = 0
                for z in range(4):
                    Elem4.suma_H[x, y] += H[z, x, y]

        print("\nMacierz H po sumowaniu:\n", Elem4.suma_H)


element1 = Elem4()
element1.pochodne()
element1.jakobian()
element1.pochodne2()
element1.macierzH()