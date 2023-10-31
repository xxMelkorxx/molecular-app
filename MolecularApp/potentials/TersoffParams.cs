using System;

namespace MolecularApp.potentials;

public struct TersoffParams
{
    /// <summary>
    /// 1 эВ в Дж.
    /// </summary>
    public const double Ev = 1.602176634e-19;

    public double A, B, S, R, b, c, d, n, h, l1, l2;

    /// <summary>
    /// Параметры потенциала для Si (кремний).
    /// </summary>
    public static TersoffParams ParamsSi => new(1830.8, 471.18, 24.799, 17.322, 1.1e-6, 0.78734, 1.0039e+5, 16.217, -0.59825, 0.27, 0.3);
    
    /// <summary>
    /// Параметры потенциала для Ge (германий).
    /// </summary>
    public static TersoffParams ParamsGe => new(1769, 419.23, 24.451, 17.047, 9.0166e-7, 0.75627, 1.0643e+5, 15.652, -0.43884, 0.28, 0.31);

    /// <summary>
    /// Параметры потенциала для Sn (олово).
    /// </summary>
    public static TersoffParams ParamsSn => new(520.4677, 281.4117, 15.5, 12.5649, 6.01e-7, 0.74, 1.4e+5, 14.5, -0.502, 0.3, 0.34);
    // public static TersoffPotential.PotentialParams ParamsSn => new TersoffPotential.PotentialParams(526.46, 296.83, 0.34, 0.30, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 15.3, 12.56);
    // public static TersoffPotential.PotentialParams ParamsSn => new TersoffPotential.PotentialParams(2848, 658.62, 0.32, 0.28, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 22.5, 16.2);
    
    /// <summary>
    /// Инициализация параметров потенциала Терсоффа. 
    /// </summary>
    /// <param name="A">(эВ)</param>
    /// <param name="B">(эВ)</param>
    /// <param name="l1">(нм^-1)</param>
    /// <param name="l2">(нм^-1)</param>
    /// <param name="b">(нм^-1)</param>
    /// <param name="n"></param>
    /// <param name="c"></param>
    /// <param name="d"></param>
    /// <param name="h"></param>
    /// <param name="R">(нм)</param>
    /// <param name="S">(нм)</param>
    public TersoffParams(double A, double B, double l1, double l2, double b, double n, double c, double d, double h, double R, double S)
    {
        this.A = A * Ev;
        this.B = B * Ev;
        this.l1 = l1 * 1e9;
        this.l2 = l2 * 1e9;
        this.b = b;
        this.n = n;
        this.c = c;
        this.d = d;
        this.h = h;
        this.R = R * 1e-9;
        this.S = S * 1e-9;
    }

    /// <summary>
    /// Инициализация параметров потенциала Терсоффа для сплава.
    /// </summary>
    /// <param name="p1"></param>
    /// <param name="p2"></param>
    public TersoffParams(TersoffParams p1, TersoffParams p2)
    {
        A = Math.Sqrt(p1.A * p2.A);
        B = Math.Sqrt(p1.B * p2.B);
        l1 = 0.5 * (p1.l1 + p2.l1);
        l2 = 0.5 * (p1.l2 + p2.l2);
        b = 0.5 * (p1.b + p2.b);
        c = 0.5 * (p1.c + p2.c);
        d = 0.5 * (p1.d + p2.d);
        n = 0.5 * (p1.n + p2.n);
        h = 0.5 * (p1.h + p2.h);
        R = Math.Sqrt(p1.R * p2.R);
        S = Math.Sqrt(p1.S * p2.S);
    }
}