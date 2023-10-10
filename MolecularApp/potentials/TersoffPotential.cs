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

    private PotentialParams _paramsSn, _paramsGe, _paramsSnGe;

    public TersoffPotential(PotentialParams paramsGe, PotentialParams paramsSn)
    {
        this._paramsGe = paramsGe;
        this._paramsSn = paramsSn;
        this._paramsSnGe = new PotentialParams(paramsGe, paramsSn);
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

    private double DzetaIJ(PotentialParams p, Atom atomI, Dictionary<PairIndexes, double> distancesIJ, double rij, int idxJ)
    {
        var dzeta = 0d;
        var atomJ = atomI.Neighbours[idxJ];
        for (var k = 0; k < atomI.Neighbours.Count; k++)
        {
            if (k == idxJ) continue;
            var atomK = atomI.Neighbours[k];

            var rik = distancesIJ[PairIndexes.GetIndexes(atomI, atomK)];
            var pik = (atomI.Type == atomK.Type) ? (atomI.Type == AtomType.Ge ? _paramsGe : _paramsSn) : _paramsSnGe;

            if (rik > pik.S)
                return 0;
            
            var rjk = distancesIJ[PairIndexes.GetIndexes(atomJ, atomK)];

            var c2 = p.c * p.c;
            var d2 = p.d * p.d;
            var cos = 1.0 + c2 / d2 - c2 / (d2 + double.Pow(p.h - Cos(rij, rik, rjk), 2));

            dzeta += Fc(pik, rik) * cos;
        }

        return dzeta;
    }

    /// <summary>
    /// Получение радиуса обрезания.
    /// </summary>
    /// <param name="fraction"></param>
    /// <returns></returns>
    public double GetRadiusCutoff(double fraction) => fraction >= 0.5 ? _paramsGe.S : _paramsSn.S;

    public double PotentialEnergy(Atom atomI, Dictionary<PairIndexes, double> distancesIJ)
    {
        var energy = 0d;
        for (var j = 0; j < atomI.Neighbours.Count; j++)
        {
            var atomJ = atomI.Neighbours[j];

            var rij = distancesIJ[PairIndexes.GetIndexes(atomI, atomJ)];
            var pij = (atomI.Type == atomJ.Type) ? (atomI.Type == AtomType.Ge ? _paramsGe : _paramsSn) : _paramsSnGe;
            var pi = atomI.Type == AtomType.Ge ? _paramsGe : _paramsSn;

            if (rij > pij.S)
                return 0;
            
            var dzeta = DzetaIJ(pi, atomI, distancesIJ, rij, j);
            var bij = double.Pow(1 + double.Pow(pi.b * dzeta, pi.n), -1d / (2 * pi.n));

            energy += Fc(pij, rij) * (Fa(pij, rij) + bij * Fr(pij, rij));
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