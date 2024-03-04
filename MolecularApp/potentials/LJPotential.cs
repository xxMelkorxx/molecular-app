namespace MolecularApp.potentials;

public class LJPotential : IPotential
{
    public double GetRadiusCutoff(double fraction)
    {
        return 0d;
    }

    /// <summary>
    /// /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <returns></returns>
    public double PotentialDerivative(Atom selAtom)
    {
        return 0d;
    }

    /// <summary>
    /// Потенциальная энергия двух атомов (Дж).
    /// </summary>
    /// <param name="selAtom"></param>
    /// <returns></returns>
    public double PotentialEnergy(Atom selAtom)
    {
        return 0d;
    }
}