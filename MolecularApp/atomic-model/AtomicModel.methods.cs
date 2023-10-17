using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MolecularApp.atomic_model;

public partial class AtomicModel
{
    /// <summary>
    /// Объект-заглушка для синхронизации потоков.
    /// </summary>
    private object _taskLocker = new();

    /// <summary>
    /// Поиск соседей для каждого атома системы/>.
    /// </summary>
    /// <param name="searchRadius">Радиус поиска (в нм).</param>
    public void SearchAtomsNeighbours(double searchRadius)
    {
        var searchRadiusSquared = searchRadius * searchRadius;
        Atoms.ForEach(atom => atom.Neighbours.Clear());
        DistanceBetweenAtoms.Clear();
        PairIndexes indexes;

        // Расчёт расстояний между атомами системы.
        Parallel.For(0, CountAtoms, i =>
        {
            var atomI = Atoms[i];
            for (var j = i + 1; j < CountAtoms; j++)
            {
                var atomJ = Atoms[j];
                var distanceSquared = SeparationSqured(atomI.Position, atomJ.Position, out _);
                if (distanceSquared >= searchRadiusSquared)
                    continue;

                var distance = double.Sqrt(distanceSquared);
                indexes = PairIndexes.GetIndexes(atomI, atomJ);

                lock (_taskLocker)
                {
                    DistanceBetweenAtoms[indexes] = distance;
                    atomI.Neighbours.Add(atomJ);
                    atomJ.Neighbours.Add(atomI);
                }
            }
        });

        // Расчёт расстояний между соседями каждого атома системы.
        Parallel.For(0, CountAtoms, n =>
        {
            var selAtom = Atoms[n];
            for (var i = 0; i < selAtom.Neighbours.Count; i++)
            for (var j = i + 1; j < selAtom.Neighbours.Count; j++)
            {
                var atomI = selAtom.Neighbours[i];
                var atomJ = selAtom.Neighbours[j];
                var distance = Separation(atomI.Position, atomJ.Position, out _);
                indexes = PairIndexes.GetIndexes(atomI, atomJ);

                lock (_taskLocker)
                    DistanceBetweenAtoms[indexes] = distance;
            }
        });
    }

    /// <summary>
    /// Начальное смещение атомов.
    /// </summary>
    /// <param name="k">Коэффициент смещения.</param>
    public void AtomsDisplacement(double k)
    {
        Atoms.ForEach(atom =>
        {
            var displacement = (-1 * XYZ.One + 2 * new XYZ(_rnd.NextDouble(), _rnd.NextDouble(), _rnd.NextDouble())) * k * LatticeGeSn;
            Flux = XYZ.Zero;
            atom.Position = Periodic(atom.Position + displacement);
            atom.PositionNp += displacement;
        });
    }

    /// <summary>
    /// Начальная перенормировка скоростей.
    /// </summary>
    /// <param name="temp">Заданная температура.</param>
    public void InitVelocityNormalization(double temp)
    {
        const double pi2 = 2 * double.Pi;
        Atoms.ForEach(atom =>
        {
            var r1 = _rnd.NextDouble();
            var r2 = _rnd.NextDouble();
            atom.Velocity = new XYZ(
                double.Sin(pi2 * r1) * double.Cos(pi2 * r2),
                double.Sin(pi2 * r1) * double.Sin(pi2 * r2),
                double.Sin(pi2 * r1)) * double.Sqrt(3 * kB * temp / atom.Weight);
        });
    }

    /// <summary>
    /// Перенормировка скоростей к заданной температуре.
    /// </summary>
    /// <param name="temp">Заданная температура</param>
    public void VelocityNormalization(double temp)
    {
        var sum = Atoms.Sum(atom => atom.Weight * atom.Velocity.SquaredMagnitude());
        if (sum == 0)
            throw new DivideByZeroException();
        var beta = Math.Sqrt(3 * CountAtoms * kB * temp / sum);
        Atoms.ForEach(atom => atom.Velocity *= beta);
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

    /// <summary>
    /// Получение радиального распределения атомов g(r(нм)).
    /// </summary>
    /// <returns>Массив точек функции радиального распределения</returns>
    public PointD[] GetRadialDistribution()
    {
        var dr = 0.05 * LatticeGeSn * 0.726;
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
            var coef = V / (CountAtoms * 4 * Math.PI * Math.PI * rd[i].X * rd[i].X * dr);
            rd[i].Y /= countAtoms == 0 ? 1 : countAtoms;
            rd[i].Y *= 1 / coef == 0 ? 1 : coef;
        }

        return rd;
    }

    /// <summary>
    /// Получение координат атомов на текущем шаге.
    /// </summary>
    public List<XYZ> GetPosAtoms() => Atoms.Select(atom => atom.Position).ToList();

    /// <summary>
    /// Получение координат атомов без учёта ПГУ на текущем шаге.
    /// </summary>
    private List<XYZ> GetPosNpAtoms() => Atoms.Select(atom => atom.PositionNp).ToList();

    /// <summary>
    /// Получение скоростей атомов.
    /// </summary>
    /// <returns></returns>
    private List<XYZ> GetVelocitiesAtoms() => Atoms.Select(atom => atom.Velocity).ToList();

    /// <summary>
    /// Вычисление среднего квадрата смещения на текущем шаге.
    /// </summary>
    /// <returns>Средний квадрат смещения (м²)</returns>
    public double GetMsd() => _rt0.Zip(GetPosNpAtoms(), (vec1, vec2) => (vec2 - vec1).SquaredMagnitude()).Sum() / CountAtoms;

    /// <summary>
    /// Расчет коэффициента самодиффузии из среднего квадрата смещения.
    /// </summary>
    /// <param name="msdPoints">Список точек среднего квадрата смещения.</param>
    /// <param name="errorRate">Погрешность коэффициента самодиффузии.</param>
    /// <returns>Коэффициент самодиффузии (м²/с).</returns>
    public double GetSelfDiffCoefFromMsd(IEnumerable<PointD> msdPoints, out double errorRate)
    {
        var msd = msdPoints.Skip(1).ToList();
        var n = msd.Count;

        // Метод наименьших квадратов.
        var sumx = msd.Sum(p => p.X);
        var sumy = msd.Sum(p => p.Y);
        var sumxy = msd.Sum(p => p.X * p.Y);
        var sumxx = msd.Sum(p => p.X * p.X);
        var averX = sumx / n;
        var b = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
        var g = (sumxx * sumy - sumx * sumxy) / (n * sumxx - sumx * sumx);
        var q = msd.Sum(p => double.Pow(b * p.X + g - p.Y, 2));
        var sb = double.Sqrt(q / ((n - 2) * msd.Sum(p => double.Pow(p.X - averX, 2))));
        errorRate = 1.96 * sb / 6;

        return b / 6;
    }

    /// <summary>
    /// Расчет коэффициента самодиффузии из среднего квадрата смещения.
    /// </summary>
    /// <param name="p1"></param>
    /// <param name="p2"></param>
    /// <returns>Коэффициент самодиффузии (м²/с)</returns>
    public double GetSelfDiffCoefFromMsd(PointD p1, PointD p2) => (p2.Y - p1.Y) / (p2.X - p1.X) / 6d;

    /// <summary>
    /// Рассчёт автокорреляционной функции скорости атомов.
    /// </summary>
    /// <returns>Массив значений АКФ скорости.</returns>
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

    public double GetSigma() => _potential.GetRadiusCutoff(FractionGe);

    /// <summary>
    /// Вычисление параметра решётки системы по закону Вегарда.
    /// </summary>
    /// <param name="fraction">Доля первого типа атома.</param>
    /// <param name="lattice1">Параметр решётки первого типа атома.</param>
    /// <param name="lattice2">Параметр решётки второго типа атома.</param>
    /// <returns></returns>
    public static double VegardLaw(double fraction, double lattice1, double lattice2) => lattice1 * fraction + lattice2 * (1 - fraction);
}