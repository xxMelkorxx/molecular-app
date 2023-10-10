﻿namespace MolecularApp.potentials;

public interface IPotential
{
    /// <summary>
    /// 1 эВ в Дж.
    /// </summary>
    public const double Ev = 1.602176634e-19;

    public double GetRadiusCutoff(double fraction);
    
    /// <summary>
    /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    /// <param name="args"></param>
    /// <returns></returns>
    public double PotentialDerivative(object[] args);

    /// <summary>
    /// Потенциальная энергия двух атомов (Дж).
    /// </summary>
    /// <param name="args"></param>
    /// <returns></returns>
    public double PotentialEnergy(object[] args);
}