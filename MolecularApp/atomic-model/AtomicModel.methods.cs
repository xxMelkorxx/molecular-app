using System;
using System.Collections.Generic;
using System.Linq;

namespace MolecularApp.atomic_model;

public abstract partial class AtomicModel
{
    /// <summary>
    /// Учёт периодических граничных условий.
    /// </summary>
    protected XYZ Periodic(XYZ pos)
    {
        XYZ newPos;
        // Ось X.               
        if (pos.X > BoxSize)
            newPos.X = pos.X - BoxSize;
        else if (pos.X < 0)
            newPos.X = pos.X + BoxSize;
        else newPos.X = pos.X;
        // Ось Y.
        if (pos.Y > BoxSize)
            newPos.Y = pos.Y - BoxSize;
        else if (pos.Y < 0)
            newPos.Y = pos.Y + BoxSize;
        else newPos.Y = pos.Y;
        // Ось Z.
        if (pos.Z > BoxSize)
            newPos.Z = pos.Z - BoxSize;
        else if (pos.Z < 0)
            newPos.Z = pos.Z + BoxSize;
        else newPos.Z = pos.Z;

        return newPos;
    }
    
    /// <summary>
    /// Учёт периодических граничных условий.
    /// </summary>
    protected XYZ Periodic(XYZ pos, XYZ p)
    {
        XYZ newPos;
        // Ось X.               
        if (pos.X > BoxSize)
        {
            newPos.X = pos.X - BoxSize;
            Flux.X += p.X;
        }
        else if (pos.X < 0)
        {
            newPos.X = pos.X + BoxSize;
            Flux.X -= p.X;
        }
        else newPos.X = pos.X;
    
        // Ось Y.
        if (pos.Y > BoxSize)
        {
            newPos.Y = pos.Y - BoxSize;
            Flux.Y += p.Y;
        }
        else if (pos.Y < 0)
        {
            newPos.Y = pos.Y + BoxSize;
            Flux.Y -= p.Y;
        }
        else newPos.Y = pos.Y;
    
        // Ось Z.
        if (pos.Z > BoxSize)
        {
            newPos.Z = pos.Z - BoxSize;
            Flux.Z += p.Z;
        }
        else if (pos.Z < 0)
        {
            newPos.Z = pos.Z + BoxSize;
            Flux.Z -= p.Z;
        }
        else newPos.Z = pos.Z;
    
        return newPos;
    }
    
    /// <summary>
    /// Вычисление расстояния между частицами с учётом периодических граничных условий.
    /// </summary>
    protected double Separation(XYZ vec1, XYZ vec2) => Math.Sqrt(SeparationSqured(vec1, vec2));

    /// <summary>
    /// Вычисление квадрата расстояния между частицами с учётом периодических граничных условий.
    /// </summary>
    protected double SeparationSqured(XYZ vec1, XYZ vec2)
    {
        var dxdydz = vec1 - vec2;

        // Обеспечивает, что расстояние между частицами никогда не будет больше L/2.
        if (Math.Abs(dxdydz.X) > 0.5 * BoxSize)
            dxdydz.X -= Math.Sign(dxdydz.X) * BoxSize;
        if (Math.Abs(dxdydz.Y) > 0.5 * BoxSize)
            dxdydz.Y -= Math.Sign(dxdydz.Y) * BoxSize;
        if (Math.Abs(dxdydz.Z) > 0.5 * BoxSize)
            dxdydz.Z -= Math.Sign(dxdydz.Z) * BoxSize;

        return dxdydz.SquaredMagnitude();
    }
    
    /// <summary>
    /// Начальное смещение атомов.
    /// </summary>
    /// <param name="k">Коэффициент смещения.</param>
    public void AtomsDisplacement(double k)
    {
        Atoms.ForEach(atom =>
        {
            var displacement = (-1 * XYZ.One + 2 * new XYZ(_rnd.NextDouble(), _rnd.NextDouble(), _rnd.NextDouble())) * k * SystemLattice;
            atom.Position = Periodic(atom.Position + displacement);
            atom.PositionNp += displacement;
        });
    }
    
    /// <summary>
    /// Получение координат атомов на текущем шаге.
    /// </summary>
    public List<AtomItem> GetAtomItems() => Atoms.Select(atom => new AtomItem(atom.Type, atom.Position)).ToList();
    
    /// <summary>
    /// Получение скоростей атомов.
    /// </summary>
    protected List<XYZ> GetVelocitiesAtoms() => Atoms.Select(atom => atom.Velocity).ToList();
    
    /// <summary>
    /// Расчет коэффициента самодиффузии из среднего квадрата смещения.
    /// </summary>
    /// <param name="msdPoints">Список точек среднего квадрата смещения.</param>
    /// <param name="errorRate">Погрешность коэффициента самодиффузии.</param>
    /// <returns>Коэффициент самодиффузии (м²/с).</returns>
    public static double GetSelfDiffCoefFromMsd(List<PointD> msdPoints, out double errorRate)
    {
        var n = msdPoints.Count;

        // Метод наименьших квадратов.
        var sumx = msdPoints.Sum(p => p.X);
        var sumy = msdPoints.Sum(p => p.Y);
        var sumxy = msdPoints.Sum(p => p.X * p.Y);
        var sumxx = msdPoints.Sum(p => p.X * p.X);
        var averX = sumx / n;
        var b = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
        var g = (sumxx * sumy - sumx * sumxy) / (n * sumxx - sumx * sumx);
        var q = msdPoints.Sum(p => Math.Pow(b * p.X + g - p.Y, 2));
        var sb = Math.Sqrt(q / ((n - 2) * msdPoints.Sum(p => Math.Pow(p.X - averX, 2))));
        errorRate = 1.96 * sb / 6d;

        return b / 6d;
    }

    /// <summary>
    /// Расчет коэффициента самодиффузии из среднего квадрата смещения.
    /// </summary>
    public static double GetSelfDiffCoefFromMsd(PointD p1, PointD p2) => (p2.Y - p1.Y) / (p2.X - p1.X) / 6d;

    /// <summary>
    /// Рассчёт автокорреляционной функции скорости атомов.
    /// </summary>
    public double[] GetAcfs(out double norm)
    {
        var zt = new double[CountNumberAcf];
        for (var i = 0; i < CountRepeatAcf; i++)
            for (var j = 0; j < CountNumberAcf; j++)
            {
                for (var k = 0; k < CountAtoms; k++)
                    zt[j] += k != 0
                        ? _vtList[i * StepRepeatAcf][k].X * _vtList[j + i * StepRepeatAcf][k].X +
                          _vtList[i * StepRepeatAcf][k].Y * _vtList[j + i * StepRepeatAcf][k].Y +
                          _vtList[i * StepRepeatAcf][k].Z * _vtList[j + i * StepRepeatAcf][k].Z
                        : _vtList[i * StepRepeatAcf][k].Magnitude();
                zt[j] /= CountAtoms;
            }

        norm = zt.Max();
        for (var i = 0; i < zt.Length; i++)
            zt[i] /= norm;

        return zt;
    }

    /// <summary>
    /// Расчет коэффициента самодиффузии из АКФ скорости.
    /// </summary>
    /// <param name="zt">АКФ скорости.</param>
    /// <param name="norm">Коэффициент нормировки</param>
    /// <returns>Коэффициент самодиффузии (м²/с)</returns>
    public double GetSelfDiffCoefFromAcf(double[] zt, double norm) => (zt.Sum() - (zt[0] + zt[zt.Length - 1]) / 2) * dt * norm / 3;
    
    /// <summary>
    /// Получение радиального распределения атомов g(r(нм)).
    /// </summary>
    public PointD[] GetRadialDistribution()
    {
        var dr = 0.05 * SystemLattice * 0.726;
        var dr2 = dr * dr;
        var rd = new PointD[(int)(BoxSize / dr)];
        for (var i = 0; i < rd.Length; i++)
            rd[i] = new PointD(i * dr, 0);

        // Подсчёт числа атомов в центральной части расчётной ячейки.
        var countAtoms = Atoms.Count(atom =>
            atom.Position.X > 0.25 * BoxSize && atom.Position.X < 0.75 * BoxSize &&
            atom.Position.Y > 0.25 * BoxSize && atom.Position.Y < 0.75 * BoxSize &&
            atom.Position.Z > 0.25 * BoxSize && atom.Position.Z < 0.75 * BoxSize);

        // Подсчёт n(r).
        foreach (var atomI in Atoms)
        foreach (var atomJ in Atoms)
        {
            if (atomJ.Equals(atomI)) continue;
            var r2 = (atomI.Position - atomJ.Position).SquaredMagnitude();
            for (var k = 0; k < rd.Length; k++)
                if (r2 > k * k * dr2 && r2 < (k + 1) * (k + 1) * dr2)
                    rd[k].Y++;
        }

        // Усреднение.
        for (var i = 0; i < rd.Length; i++)
        {
            var coef = GetVolume / (CountAtoms * 4 * Math.PI * Math.PI * rd[i].X * rd[i].X * dr);
            rd[i].Y /= countAtoms == 0 ? 1 : countAtoms;
            rd[i].Y *= 1 / coef == 0 ? 1 : coef;
        }

        return rd;
    }
    
    /// <summary>
    /// Начальная перенормировка скоростей.
    /// </summary>
    /// <param name="temp">Заданная температура.</param>
    public void InitVelocityNormalization(double temp)
    {
        const double pi2 = 2 * Math.PI;
        Atoms.ForEach(atom =>
        {
            var r1 = _rnd.NextDouble();
            var r2 = _rnd.NextDouble();
            atom.Velocity = new XYZ(
                Math.Sin(pi2 * r1) * Math.Cos(pi2 * r2),
                Math.Sin(pi2 * r1) * Math.Sin(pi2 * r2),
                Math.Sin(pi2 * r1)
            ) * Math.Sqrt(3 * kB * temp / atom.Weight);
        });
        PulseZeroing();
    }

    /// <summary>
    /// Перенормировка скоростей к заданной температуре.
    /// </summary>
    /// <param name="temp">Заданная температура</param>
    public void VelocityNormalization(double temp)
    {
        var sumKE = Atoms.Sum(atom => atom.Weight * atom.Velocity.SquaredMagnitude());
        if (sumKE == 0)
            throw new DivideByZeroException();

        var beta = Math.Sqrt(3 * CountAtoms * kB * temp / sumKE);
        Atoms.ForEach(atom => atom.Velocity *= beta);
        PulseZeroing();
    }
    
    /// <summary>
    /// Зануление импульса системы.
    /// </summary>
    /// <param name="eps">Точность.</param>
    public void PulseZeroing(double eps = 1e-5)
    {
        XYZ sum;
        while (true)
        {
            sum = XYZ.Zero;
            Atoms.ForEach(atom => sum += atom.Velocity);
            sum /= CountAtoms;

            if (Math.Abs(sum.X + sum.Y + sum.Z) > eps)
                Atoms.ForEach(atom => atom.Velocity -= sum);
            else break;
        }
    }
}