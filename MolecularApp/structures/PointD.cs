﻿using System;

namespace MolecularApp;

public struct PointD
{
    public double X, Y;

    public PointD(double x, double y)
    {
        X = x;
        Y = y;
    }

    public PointD Rotate(double angle) => new(X * Math.Cos(angle) - Y * Math.Sin(angle), X * Math.Sin(angle) + Y * Math.Cos(angle));
}