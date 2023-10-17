using System;
using System.Collections.Generic;
using System.Linq;

namespace MolecularApp.potentials;

public class TersoffPotential : IPotential
{
    private TersoffParams _paramsSn, _paramsGe, _paramsSnGe;

    /// <summary>
    /// Инициализация потенциала.
    /// </summary>
    /// <param name="paramsGe"></param>
    /// <param name="paramsSn"></param>
    public TersoffPotential(TersoffParams paramsGe, TersoffParams paramsSn)
    {
        _paramsGe = paramsGe;
        _paramsSn = paramsSn;
        _paramsSnGe = new TersoffParams(paramsGe, paramsSn);
    }

    /// <summary>
    /// Получение радиуса обрезания.
    /// </summary>
    /// <param name="fraction"></param>
    /// <returns></returns>
    public double GetRadiusCutoff(double fraction) => fraction >= 0.5 ? _paramsGe.S : _paramsSn.S;

    /// <summary>
    /// Потенциальная энергия (Дж).
    /// </summary>
    /// <param name="atomI"></param>
    /// <param name="distancesIJ"></param>
    /// <returns></returns>
    public double PotentialEnergy(Atom atomI, Dictionary<PairIndexes, double> distancesIJ)
    {
        var energy = 0d;
        for (var j = 0; j < atomI.Neighbours.Count; j++)
        {
            var atomJ = atomI.Neighbours[j];

            var pij = (atomI.Type == atomJ.Type) ? (atomI.Type == AtomType.Ge ? _paramsGe : _paramsSn) : _paramsSnGe;
            var pi = atomI.Type == AtomType.Ge ? _paramsGe : _paramsSn;

            var rij = distancesIJ[PairIndexes.GetIndexes(atomI, atomJ)];
            if (rij > pij.S)
                return 0;

            var dzeta = DzetaIJ(pi, atomI, distancesIJ, rij, j);
            var bij = Math.Pow(1 + Math.Pow(pi.b * dzeta, pi.n), -1d / (2 * pi.n));

            energy += Fc(pij, rij) * (Fa(pij, rij) + bij * Fr(pij, rij));
        }

        return energy / 2;
    }

    /// <summary>
    /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <param name="distancesIJ"></param>
    /// <returns></returns>
    public double PotentialDerivative(Atom selAtom, Dictionary<PairIndexes, double> distancesIJ) => 
        PotentialEnergy(selAtom, distancesIJ) + selAtom.Neighbours.Sum(atom => PotentialEnergy(atom, distancesIJ));

    /// <summary>
    /// Функция обрезания.
    /// </summary>
    /// <param name="p">Параметры потенциала</param>
    /// <param name="rij">Расстояние между атомами</param>
    /// <returns></returns>
    private static double Fc(TersoffParams p, double rij) => (rij <= p.R) ? 1 : (rij >= p.S) ? 0 : 0.5 + 0.5 * Math.Cos(Math.PI * (rij - p.R) / (p.S - p.R));

    /// <summary>
    /// Потенциал притяжения.
    /// </summary>
    /// <param name="p">Параметры потенциала</param>
    /// <param name="rij">Расстояние между атомами</param>
    /// <returns></returns>
    private static double Fa(TersoffParams p, double rij) => p.A * Math.Exp(-p.l1 * rij);

    /// <summary>
    /// Потенциал отталкивания.
    /// </summary>
    /// <param name="p">Параметры потенциала</param>
    /// <param name="rij">Расстояние между атомами</param>
    /// <returns></returns>
    private static double Fr(TersoffParams p, double rij) => -p.B * Math.Exp(-p.l2 * rij);

    /// <summary>
    /// Косинус угла между связями i-j и i-k, вычисляемый по теореме косинусов.
    /// </summary>
    /// <param name="rij">Расстояние между i-j атомами</param>
    /// <param name="rik">Расстояние между i-k атомами</param>
    /// <param name="rjk">Расстояние между j-k атомами</param>
    /// <returns></returns>
    private static double Cos(double rij, double rik, double rjk) => (rik * rik + rij * rij - rjk * rjk) / (2d * rij * rik);

    private double DzetaIJ(TersoffParams p, Atom atomI, Dictionary<PairIndexes, double> distancesIJ, double rij, int idxJ)
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
            var cos = 1.0 + c2 / d2 - c2 / (d2 + Math.Pow(p.h - Cos(rij, rik, rjk), 2));

            dzeta += Fc(pik, rik) * cos;
        }

        return dzeta;
    }
}