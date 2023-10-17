using System;

namespace MolecularApp.potentials;

public struct TersoffParams
{
    public double A, B, S, R, b, c, d, n, h, l1, l2;

    /// <summary>
    /// Инициализация параметров потенциала Терсоффа. 
    /// </summary>
    /// <param name="A">эВ</param>
    /// <param name="B">эВ</param>
    /// <param name="S"></param>
    /// <param name="R"></param>
    /// <param name="b"></param>
    /// <param name="c"></param>
    /// <param name="d"></param>
    /// <param name="n"></param>
    /// <param name="h"></param>
    /// <param name="l1"></param>
    /// <param name="l2"></param>
    public TersoffParams(double A, double B, double S, double R, double b, double c, double d, double n, double h, double l1, double l2)
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

    /// <summary>
    /// Инициализация параметров потенциала Терсоффа для сплава.
    /// </summary>
    /// <param name="p1"></param>
    /// <param name="p2"></param>
    public TersoffParams(TersoffParams p1, TersoffParams p2)
    {
        A = Math.Sqrt(p1.A * p2.A);
        B = Math.Sqrt(p1.B * p2.B);
        S = Math.Sqrt(p1.S * p2.S);
        R = Math.Sqrt(p1.R * p2.R);
        b = 0.5 * (p1.b + p2.b);
        c = 0.5 * (p1.c + p2.c);
        d = 0.5 * (p1.d + p2.d);
        n = 0.5 * (p1.n + p2.n);
        h = 0.5 * (p1.h + p2.h);
        l1 = 0.5 * (p1.l1 + p2.l1);
        l2 = 0.5 * (p1.l2 + p2.l2);
    }
}