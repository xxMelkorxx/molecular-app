using System;

namespace MolecularApp.potentials;

public class MLJPotential
{
    // 1 эВ в Дж.
    public const double Ev = 1.602176634e-19;
    
    // Тип атома.
    public AtomType AtomType
    {
        set
        {
            _atomType = value;
            switch (_atomType)
            {
                case AtomType.Ar:
                    D = 0.01029 * Ev;
                    Sigma = 0.3408e-9;
                    break;
                default: throw new Exception("Неверный тип атома");
            }
        }
    }
    private AtomType _atomType;
    
    // Модуль потенциальной энергии взаимодействия между атомами при равновесии (м).
    public double Sigma;

    // Равновесное расстояние между центрами атомов
    public double R0 => Sigma * Math.Pow(2, 1d / 6d);

    // Ближний радиус обрезания потенциала.
    public double R1 => 1.2 * R0;

    // Дальний радиус обрезания потенциала.
    public double R2 => 1.8 * R0;

    // Модуль потенциальной энергии взаимодействия между атомами при равновесии (Дж).
    public double D;

    public MLJPotential(AtomType atomType)
    {
        _atomType = atomType;
    }
    
    /// <summary>
    /// Межатомная сила взаимодействия в потенциале (Дж * м).
    /// </summary>
    public XYZ PotentialDerivative(double r, XYZ dxdydz)
    {
        return (r < R1) ? Flj(r) * dxdydz : (r > R2) ? XYZ.Zero : Flj(r) * dxdydz * K(r);
    }

    /// <summary>
    /// Потенциальная энергия двух атомов (Дж).
    /// </summary>
    public double PotentialEnergy(double r2)
    {
        return r2 < R1 ? Plj(r2) : r2 > R2 ? 0 : Plj(r2) * K(double.Sqrt(r2));;
    }
    
    /// <summary>
    /// Функция обрезания потенциала.
    /// </summary>
    /// <param name="r">Расстояние между частицами.</param>
    /// <returns></returns>
    private double K(double r) => Math.Pow(1 - (r - R1) * (r - R1) / (R1 - R2) / (R1 - R2), 2);

    /// <summary>
    /// Потенциал Леннарда-Джонса.
    /// </summary>
    /// <param name="r2">Расстояние между частицами.</param>
    /// <returns></returns>
    private double Plj(double r2)
    {
        if (r2 == 0)
            throw new DivideByZeroException();
        
        var ri2 = Sigma * Sigma / r2;
        var ri6 = ri2 * ri2 * ri2;
        
        return 4 * D * ri6 * (ri6 - 1);
    }

    /// <summary>
    /// Cила в потенциале Леннарда-Джонса.
    /// </summary>
    /// <param name="r2">Расстояние между частицами.</param>
    /// <returns></returns>
    private double Flj(double r2)
    {
        if (r2 == 0)
            throw new DivideByZeroException();
        
        var ri2 = Sigma * Sigma / r2;
        var ri6 = ri2 * ri2 * ri2;

        return 24 * D * ri6 * (2 * ri6 - 1) / r2;
    }
}