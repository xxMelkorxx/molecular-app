using System.Collections.Generic;

namespace MolecularApp.potentials;

public interface IPotential
{
    public double GetRadiusCutoff(double fraction);
    
    /// <summary>
    /// /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <param name="atomsDistances"></param>
    /// <returns></returns>
    public double PotentialDerivative(Atom selAtom, Dictionary<PairIndexes, double> atomsDistances);

    /// <summary>
    /// Потенциальная энергия двух атомов (Дж).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <param name="atomsDistances"></param>
    /// <returns></returns>
    public double PotentialEnergy(Atom selAtom, Dictionary<PairIndexes, double> atomsDistances);
}