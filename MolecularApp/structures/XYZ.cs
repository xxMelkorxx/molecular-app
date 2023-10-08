using System;

namespace MolecularApp;

public struct XYZ
{
    public double X, Y, Z;

    public XYZ(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    /// <summary>
    /// Нулевой вектор.
    /// </summary>
    public static XYZ Zero => new(0, 0, 0);

    /// <summary>
    /// Единичный вектор.
    /// </summary>
    public static XYZ One => new(1, 1, 1);

    /// <summary>
    /// Возвращает значение наибольшей координаты.
    /// </summary>
    public double MaxElement() => double.Max(double.Max(X, Y), double.Max(X, Z));

    /// <summary>
    /// Возвращает значение наименьшей координаты.
    /// </summary>
    public double MinElement() => double.Min(double.Min(X, Y), double.Min(X, Z));

    /// <summary>
    /// Квадрат величины вектора.
    /// </summary>
    public double SquaredMagnitude() => X * X + Y * Y + Z * Z;

    /// <summary>
    /// Величина вектора.
    /// </summary>
    public double Magnitude() => double.Sqrt(SquaredMagnitude());

    public static XYZ operator +(XYZ vec1, XYZ vec2) => new(vec1.X + vec2.X, vec1.Y + vec2.Y, vec1.Z + vec2.Z);

    public static XYZ operator +(XYZ vec1, int value) => new(vec1.X + value, vec1.Y + value, vec1.Z + value);

    public static XYZ operator +(XYZ vec1, double value) => new(vec1.X + value, vec1.Y + value, vec1.Z + value);

    public static XYZ operator -(XYZ vec1, XYZ vec2) => new(vec1.X - vec2.X, vec1.Y - vec2.Y, vec1.Z - vec2.Z);

    public static XYZ operator -(XYZ vec1, int value) => new(vec1.X - value, vec1.Y - value, vec1.Z - value);

    public static XYZ operator -(XYZ vec1, double value) => new(vec1.X - value, vec1.Y - value, vec1.Z - value);

    public static XYZ operator *(XYZ vec, int num) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static XYZ operator *(int num, XYZ vec) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static XYZ operator *(XYZ vec, double num) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static XYZ operator *(double num, XYZ vec) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static XYZ operator /(XYZ vec, int num) => num != 0 ? new XYZ(vec.X / num, vec.Y / num, vec.Z / num) : throw new DivideByZeroException();

    public static XYZ operator /(XYZ vec, double num) => num != 0 ? new XYZ(vec.X / num, vec.Y / num, vec.Z / num) : throw new DivideByZeroException();
}