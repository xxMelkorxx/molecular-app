using System;
using System.Collections.Generic;

namespace MolecularApp.potentials;

public class TersoffPotential : IPotential
{
    public struct PotentialParams
    {
        public double A, B, S, R, b, c, d, n, h, l1, l2;

        public PotentialParams(double A, double B, double S, double R, double b, double c, double d, double n, double h, double l1, double l2)
        {
            this.A = A;
            this.B = B;
            this.S = S;
            this.R = R;
            this.b = b;
            this.c = c;
            this.d = d;
            this.n = n;
            this.h = h;
            this.l1 = l1;
            this.l2 = l2;
        }

        public PotentialParams(PotentialParams p1, PotentialParams p2)
        {
            this.A = double.Sqrt(p1.A * p2.A);
            this.B = double.Sqrt(p1.B * p2.B);
            this.S = double.Sqrt(p1.S * p2.S);
            this.R = double.Sqrt(p1.R * p2.R);
            this.b = 0.5 * (p1.b + p2.b);
            this.c = 0.5 * (p1.c + p2.c);
            this.d = 0.5 * (p1.d + p2.d);
            this.n = 0.5 * (p1.n + p2.n);
            this.h = 0.5 * (p1.h + p2.h);
            this.l1 = 0.5 * (p1.l1 + p2.l1);
            this.l2 = 0.5 * (p1.l2 + p2.l2);
        }
    }

    public PotentialParams paramsSn, paramsGe, paramsSnGe;
    // private AtomType potAtomType1, potAtomType2;

    public TersoffPotential(PotentialParams paramsGe, PotentialParams paramsSn, AtomType atomType1, AtomType atomType2)
    {
        // potAtomType1 = atomType1;
        // potAtomType2 = atomType2;

        this.paramsGe = paramsGe;
        this.paramsSn = paramsSn;
        this.paramsSnGe = new PotentialParams(paramsGe, paramsSn);
    }

    /// <summary>
    /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    /// <param name="args"></param>
    /// <returns></returns>
    public double PotentialDerivative(object[] args)
    {
        return 0;
    }

    /// <summary>
    /// Потенциальная энергия (Дж).
    /// </summary>
    /// <param name="args"></param>
    /// <returns></returns>
    public double PotentialEnergy(object[] args)
    {
        return 0;
    }

    private static double Fc(PotentialParams p, double rij)
    {
        if (rij <= p.R) return 1;
        if (rij >= p.S) return 0;
        return 0.5 + 0.5 * double.Cos(double.Pi * (rij - p.R) / (p.S - p.R));
    }

    /// <summary>
    /// Потенциал притяжения.
    /// </summary>
    /// <param name="p">Параметры потенциала</param>
    /// <param name="rij">Расстояние между атомами</param>
    /// <returns></returns>
    private static double Fa(PotentialParams p, double rij) => p.A * double.Exp(-p.l1 * rij);

    /// <summary>
    /// Потенциал отталкивания.
    /// </summary>
    /// <param name="p">Параметры потенциала</param>
    /// <param name="rij">Расстояние между атомами</param>
    /// <returns></returns>
    private static double Fr(PotentialParams p, double rij) => -p.B * double.Exp(-p.l2 * rij);

    /// <summary>
    /// Косинус угла между связями i-j и i-k, вычисляемый по теореме косинусов.
    /// </summary>
    /// <param name="Rij">Расстояние между i-j атомами</param>
    /// <param name="Rik">Расстояние между i-k атомами</param>
    /// <param name="Rjk">Расстояние между j-k атомами</param>
    /// <returns></returns>
    private static double Cos(double Rij, double Rik, double Rjk) => (double.Pow(Rik, 2) + double.Pow(Rij, 2) - double.Pow(Rjk, 2)) / (2.0 * Rij * Rik);

    private double DzetaIJ(PotentialParams p, Atom selAtom, Dictionary<PairIndexes, double> distancesIJ, double rij, int idxJ)
    {
        double dzeta = 0;
        PairIndexes indexes;
        var neighJ = selAtom.Neighbours[idxJ];

        for (int k = 0; k < selAtom.Neighbours.Count; k++)
        {
            if (k == idxJ) continue;

            Atom neighK = selAtom.Neighbours[k];

            indexes = PairIndexes.GetIndexes(selAtom, neighK);
            double rik = distancesIJ[indexes];

            PotentialParams pik = potType1;
            if (selAtom.AtomType == potAtomType2) pik = potType2;
            if (selAtom.AtomType != neighK.AtomType) pik = potType12;

            if (rik < pik.S)
            {
                indexes = PairIndexes.GetIndexes(neighJ, neighK);
                double rjk = distancesIJ[indexes];

                double cos = Cos(rij, rik, rjk);
                double c2 = p.c * p.c;
                double d2 = p.d * p.d;
                cos = 1.0 + (c2 / d2) - (c2 / (d2 + Math.Pow(p.h - cos, 2)));

                dzeta += Fc(pik, rik) * cos; //* Math.Exp(Math.Pow(pij.l2 * (rij - rik), 3));
            }
        }

        return dzeta;
    }

    /// <summary>
    /// Получение радиуса обрезания.
    /// </summary>
    /// <param name="fraction"></param>
    /// <returns></returns>
    public double GetRadiusCutoff(double fraction) => fraction >= 0.5 ? paramsGe.S : paramsSn.S;

    public double PotentialEnergy(Atom selAtom, Dictionary<PairIndexes, double> distancesIJ)
    {
        double energy = 0;
        PairIndexes indexes;

        for (int j = 0; j < selAtom.Neighbours.Count; j++)
        {
            Atom neighJ = selAtom.Neighbours[j];

            indexes = PairIndexes.GetIndexes(selAtom, neighJ);
            double rij = distancesIJ[indexes];

            PotentialParams pij = potType1;
            if (selAtom.AtomType == potAtomType2) pij = potType2;
            if (selAtom.AtomType != neighJ.AtomType) pij = potType12;

            PotentialParams pi = potType1;
            if (selAtom.AtomType == potAtomType2) pi = potType2;

            if (rij < pij.S)
            {
                double dzeta = DzetaIJ(pi, selAtom, distancesIJ, rij, j);
                double bij = Math.Pow(1 + Math.Pow(pi.b * dzeta, pi.n), -1.0 / (2 * pi.n)); //-1.0 / 2 * pi.n

                energy += Fc(pij, rij) * (Fa(pij, rij) + bij * Fr(pij, rij));
            }
        }

        return energy / 2;
    }

    public double PotentialDerivative(Atom selAtom, Dictionary<PairIndexes, double> distancesIJ)
    {
        double energy = PotentialEnergy(selAtom, distancesIJ);

        for (int i = 0; i < selAtom.Neighbours.Count; i++)
        {
            Atom neigh = selAtom.Neighbours[i];
            double atomEnergy = PotentialEnergy(neigh, distancesIJ);
            energy += atomEnergy;
        }

        return energy;
    }
}