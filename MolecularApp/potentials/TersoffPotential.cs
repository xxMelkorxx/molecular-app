using System;
using System.Collections.Generic;
using System.Linq;

namespace MolecularApp.potentials;

public class TersoffPotential : IPotential
{
    private AtomType _firstTypeAtom, _secondTypeAtom;
    private TersoffParams _secondAtomParams, _firstAtomParams, _commonAtomParams;

    /// <summary>
    /// Инициализация потенциала.
    /// </summary>
    /// <param name="firstTypeAtom"></param>
    /// <param name="secondTypeAtom"></param>
    public TersoffPotential(AtomType firstTypeAtom, AtomType secondTypeAtom)
    {
        _firstTypeAtom = firstTypeAtom;
        _secondTypeAtom = secondTypeAtom;

        _firstAtomParams = firstTypeAtom switch
        {
            AtomType.Si => TersoffParams.ParamsSi,
            AtomType.Ge => TersoffParams.ParamsGe,
            AtomType.Sn => TersoffParams.ParamsSn,
            _ => throw new ArgumentException("Неверный тип атома")
        };
        _secondAtomParams = secondTypeAtom switch
        {
            AtomType.Si => TersoffParams.ParamsSi,
            AtomType.Ge => TersoffParams.ParamsGe,
            AtomType.Sn => TersoffParams.ParamsSn,
            _ => throw new ArgumentException("Неверный тип атома")
        };
        _commonAtomParams = new TersoffParams(_firstAtomParams, _secondAtomParams);
    }

    /// <summary>
    /// Получение радиуса обрезания.
    /// </summary>
    /// <param name="fraction"></param>
    /// <returns></returns>
    public double GetRadiusCutoff(double fraction) => fraction >= 0.5 ? _firstAtomParams.S : _secondAtomParams.S;

    /// <summary>
    /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <param name="atomsDistances"></param>
    /// <returns></returns>
    public double PotentialDerivative(Atom selAtom, Dictionary<PairIndexes, double> atomsDistances) =>
        PotentialEnergy(selAtom, atomsDistances) + selAtom.Neighbours.Sum(atom => PotentialEnergy(atom, atomsDistances));

    /// <summary>
    /// Потенциальная энергия (Дж).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <param name="atomsDistances"></param>
    /// <returns></returns>
    public double PotentialEnergy(Atom selAtom, Dictionary<PairIndexes, double> atomsDistances)
    {
        var energy = 0d;
        for (var j = 0; j < selAtom.Neighbours.Count; j++)
        {
            var neigh = selAtom.Neighbours[j];

            var paramsIJ = (selAtom.Type == neigh.Type) ? (selAtom.Type == _firstTypeAtom ? _firstAtomParams : _secondAtomParams) : _commonAtomParams;
            var paramsI = selAtom.Type == _firstTypeAtom ? _firstAtomParams : _secondAtomParams;
            var rij = atomsDistances[PairIndexes.GetIndexes(selAtom, neigh)];

            if (rij < paramsIJ.S)
            {
                var dzetaIJ = DzetaIJ(selAtom, neigh, paramsI, atomsDistances, rij);
                var bij = Math.Pow(1 + Math.Pow(paramsI.b * dzetaIJ, paramsI.n), -1d / (2 * paramsI.n));

                energy += Fc(paramsIJ, rij) * (Fa(paramsIJ, rij) + bij * Fr(paramsIJ, rij));
            }
        }

        return energy / 2;
    }

    private double DzetaIJ(Atom atomI, Atom atomJ, TersoffParams p, IReadOnlyDictionary<PairIndexes, double> atomsDistances, double rij)
    {
        var dzetaIJ = 0d;
        for (var k = 0; k < atomI.Neighbours.Count; k++)
        {
            var atomK = atomI.Neighbours[k];
            if (atomK.Index == atomJ.Index) continue;

            var rik = atomsDistances[PairIndexes.GetIndexes(atomI, atomK)];
            var paramsIK = (atomI.Type == atomK.Type) ? (atomI.Type == _firstTypeAtom ? _firstAtomParams : _secondAtomParams) : _commonAtomParams;

            if (rik < paramsIK.S)
            {
                var rjk = atomsDistances[PairIndexes.GetIndexes(atomJ, atomK)];
                var gijk = 1d + p.c * p.c / (p.d * p.d) - p.c * p.c / (p.d * p.d + Math.Pow(p.h - CosOijk(rij, rik, rjk), 2));
                dzetaIJ += Fc(paramsIK, rik) * gijk;
            }
        }

        return dzetaIJ;
    }

    /// <summary>
    /// Косинус угла между связями i-j и i-k, вычисляемый по теореме косинусов.
    /// </summary>
    /// <param name="rij">Расстояние между i-j атомами</param>
    /// <param name="rik">Расстояние между i-k атомами</param>
    /// <param name="rjk">Расстояние между j-k атомами</param>
    /// <returns></returns>
    private static double CosOijk(double rij, double rik, double rjk) => (rik * rik + rij * rij - rjk * rjk) / (2d * rij * rik);

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
}