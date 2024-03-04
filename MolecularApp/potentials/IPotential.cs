namespace MolecularApp.potentials;

public interface IPotential
{
    /// <summary>
    /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    public double PotentialDerivative(Atom selAtom);

    /// <summary>
    /// Потенциальная энергия двух атомов (Дж).
    /// </summary>
    public double PotentialEnergy(Atom selAtom);
}