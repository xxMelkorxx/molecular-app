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