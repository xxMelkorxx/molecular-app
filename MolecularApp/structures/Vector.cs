using System;

namespace MolecularApp;

public struct Vector
{
    public double X, Y, Z;

    public Vector(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public static Vector Zero => new(0, 0, 0);

    public static Vector One => new(1, 1, 1);

    public double MaxElement() => Math.Max(Math.Max(X, Y), Math.Max(X, Z));

    public double MinElement() => Math.Min(Math.Min(X, Y), Math.Min(X, Z));

    public double SquaredMagnitude() => X * X + Y * Y + Z * Z;

    public double Magnitude() => Math.Sqrt(SquaredMagnitude());

    public static Vector operator +(Vector vec1, Vector vec2) => new(vec1.X + vec2.X, vec1.Y + vec2.Y, vec1.Z + vec2.Z);

    public static Vector operator +(Vector vec1, int value) => new(vec1.X + value, vec1.Y + value, vec1.Z + value);

    public static Vector operator +(Vector vec1, double value) => new(vec1.X + value, vec1.Y + value, vec1.Z + value);

    public static Vector operator -(Vector vec1, Vector vec2) => new(vec1.X - vec2.X, vec1.Y - vec2.Y, vec1.Z - vec2.Z);

    public static Vector operator -(Vector vec1, int value) => new(vec1.X - value, vec1.Y - value, vec1.Z - value);

    public static Vector operator -(Vector vec1, double value) => new(vec1.X - value, vec1.Y - value, vec1.Z - value);

    public static Vector operator *(Vector vec, int num) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static Vector operator *(int num, Vector vec) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static Vector operator *(Vector vec, double num) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static Vector operator *(double num, Vector vec) => new(vec.X * num, vec.Y * num, vec.Z * num);

    public static Vector operator /(Vector vec, int num) => num != 0 ? new Vector(vec.X / num, vec.Y / num, vec.Z / num) : throw new DivideByZeroException();

    public static Vector operator /(Vector vec, double num) => num != 0 ? new Vector(vec.X / num, vec.Y / num, vec.Z / num) : throw new DivideByZeroException();
}