using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MolecularApp.atomic_model;

public partial class AtomicModel
{
    /// <summary>
    /// Вычисление начальных параметров системы (Acceleration, Pe, Ke, Press).
    /// </summary>
    public void InitCalculation()
    {
        // Рассчёт радиуса поиска соседей с учётом величины сдвига при расчёте ускорений (методом трёхточечного дифференцирования).
        var searchRadius = _potential.GetRadiusCutoff(FractionGe) + 1.2 * (1e-3 * LatticeGeSn);
        // Рассчёт соседей для каждого атома.
        SearchAtomsNeighbours(searchRadius);
        // Начальный подсчёт ускорений атомов.
        Accels();
    }

    /// <summary>
    /// Алгоритм Верле для вычисления координат и скоростей атомов на временном шаге.
    /// </summary>
    /// <param name="searchRadius">Радиус поиска (в нм).</param>
    public void Verlet(double searchRadius)
    {
        Atoms.ForEach(atom =>
        {
            var newPos = atom.Velocity * dt + 0.5 * atom.Acceleration * dt * dt;
            atom.Position = Periodic(atom.Position + newPos);
            atom.PositionNp += newPos;
        });
        // Поиск соседей атомов после их сдвига.
        SearchAtomsNeighbours(searchRadius);
        // Расчёт новых скоростей и ускорений.
        Atoms.ForEach(atom => atom.Velocity = 0.5 * atom.Acceleration * dt);
        Accels();
        Atoms.ForEach(atom => atom.Velocity = 0.5 * atom.Acceleration * dt);
    }

    /// <summary>
    /// Подсчёт ускорения системы атомов.
    /// </summary>
    public void Accels()
    {
        Parallel.For(0, CountAtoms, i =>
        {
            // Считаем, что при сдивиге у атома соседи не меняются.
            _distIJAround.Clear();
            
            var delta = 1e-3 * LatticeGeSn;
            var force = ThreePointsDiffAtom(Atoms[i], delta, ref _distIJAround);

            // эВ -> Дж -> а.е.м. * нм2/пс2 => эВ*96.352
            const double eVToDzh = 96.4830626;
            force = new XYZ(-force.X, -force.Y, -force.Z) * eVToDzh;

            Atoms[i].Acceleration = force / Atoms[i].Weight;
        });
    }

    /// <summary>
    /// Расчёт производной потенциальной энергии по формуле трехточечного дифференцирования.
    /// </summary>
    /// <param name="a">Выбранный атом.</param>
    /// <param name="delta">Сдвиг атома по координате.</param>
    /// <param name="distIJAround">Словарь, в который записываются расстояния между атомами.</param>
    /// <returns></returns>
    private XYZ ThreePointsDiffAtom(Atom a, double delta, ref Dictionary<PairIndexes, double> distIJAround)
    {
        double potentialLeft, potentialRight;
        XYZ force;

        potentialRight = PotentialEnergyAtomShifted(a, Periodic(a.Position + new XYZ(delta, 0, 0)), ref distIJAround);
        potentialLeft = PotentialEnergyAtomShifted(a, Periodic(a.Position - new XYZ(delta, 0, 0)), ref distIJAround);
        force.X = potentialRight - potentialLeft;

        potentialRight = PotentialEnergyAtomShifted(a, Periodic(a.Position + new XYZ(0, delta, 0)), ref distIJAround);
        potentialLeft = PotentialEnergyAtomShifted(a, Periodic(a.Position - new XYZ(0, delta, 0)), ref distIJAround);
        force.Y = potentialRight - potentialLeft;

        potentialRight = PotentialEnergyAtomShifted(a, Periodic(a.Position + new XYZ(0, 0, delta)), ref distIJAround);
        potentialLeft = PotentialEnergyAtomShifted(a, Periodic(a.Position - new XYZ(0, 0, delta)), ref distIJAround);
        force.Z = potentialRight - potentialLeft;

        return force / (2 * delta);
    }

    /// <summary>
    /// Подсчёт потенциальной энергии атома и его окружения при сдвиге.
    /// </summary>
    /// <param name="a">Выбранный атом.</param>
    /// <param name="shiftedPosition">Сдвинутая координата атома.</param>
    /// <param name="distIJAround">Словарь, в который записываются расстояния между атомами.</param>
    /// <returns></returns>
    private double PotentialEnergyAtomShifted(Atom a, XYZ shiftedPosition, ref Dictionary<PairIndexes, double> distIJAround)
    {
        RecalculateDistanceToNeighbours(shiftedPosition, a, ref distIJAround);
        return _potential.PotentialDerivative(a, distIJAround);
    }

    /// <summary>
    /// Подсчитывает расстояния между атомами, а именно - выбранного атома, его первый и вторых соседей;
    /// считаем что при сдвиге атома его соседи не изменяются!
    /// </summary>
    /// <param name="shiftedPosition">Сдвинутая координата атома.</param>
    /// <param name="atom">Сдвигаемый атом.</param>
    /// <param name="distanceIJAroundSelAtom">Словарь, в который записываются расстояния между атомами.</param>
    private void RecalculateDistanceToNeighbours(XYZ shiftedPosition, Atom atom, ref Dictionary<PairIndexes, double> distanceIJAroundSelAtom)
    {
        for (var n = 0; n < atom.Neighbours.Count; n++)
        {
            var neigh = atom.Neighbours[n];
            var distance = Separation(shiftedPosition, neigh.Position, out _);
            distanceIJAroundSelAtom[PairIndexes.GetIndexes(atom, neigh)] = distance;

            for (var k = 0; k < neigh.Neighbours.Count; k++)
            {
                var neighK = neigh.Neighbours[k];
                if (neighK.Index == atom.Index) continue;

                var indexes = PairIndexes.GetIndexes(neigh, neighK);
                distanceIJAroundSelAtom[indexes] = DistanceBetweenAtoms[indexes];
            }
        }

        // Вычисляем расстояние для угловых частей потенциалов.
        RecalcDistanceJK(shiftedPosition, atom, atom, ref distanceIJAroundSelAtom);
        for (var n = 0; n < atom.Neighbours.Count; n++)
        {
            var neigh = atom.Neighbours[n];
            RecalcDistanceJK(shiftedPosition, atom, neigh, ref distanceIJAroundSelAtom);
        }
    }

    /// <summary>
    /// Расчёт расстояний между соседями выбранного атома.
    /// </summary>
    /// <param name="shiftedPosition">Сдвинутая координата атома.</param>
    /// <param name="shiftedAtom">Сдвигаемый атом.</param>
    /// <param name="selAtom">Выбранный атом.</param>
    /// <param name="atomsDistanceJKAround">Словарь, в который записываются расстояния между соседями выбранного атома.</param>
    private void RecalcDistanceJK(XYZ shiftedPosition, Atom shiftedAtom, Atom selAtom, ref Dictionary<PairIndexes, double> atomsDistanceJKAround)
    {
        for (var i = 0; i < selAtom.Neighbours.Count; i++)
        for (var j = i + 1; j < selAtom.Neighbours.Count; j++)
        {
            var neighI = selAtom.Neighbours[i];
            var neighJ = selAtom.Neighbours[j];
            var indexes = PairIndexes.GetIndexes(neighI, neighJ);

            double distance;
            if (neighI.Index == shiftedAtom.Index)
                distance = Separation(shiftedPosition, neighJ.Position, out _);
            else if (neighJ.Index == shiftedAtom.Index)
                distance = Separation(neighI.Position, shiftedPosition, out _);
            else
                distance = DistanceBetweenAtoms[indexes];

            atomsDistanceJKAround[indexes] = distance;
        }
    }

    /// <summary>
    /// Вычисление расстояния между частицами с учётом периодических граничных условий. 
    /// </summary>
    /// <param name="vec1"></param>
    /// <param name="vec2"></param>
    /// <param name="dxdydz"></param>
    /// <returns></returns>
    private double Separation(XYZ vec1, XYZ vec2, out XYZ dxdydz) => Math.Sqrt(SeparationSqured(vec1, vec2, out dxdydz));

    /// <summary>
    /// Вычисление квадрата расстояния между частицами с учётом периодических граничных условий. 
    /// </summary>
    /// <param name="vec1"></param>
    /// <param name="vec2"></param>
    /// <param name="dxdydz"></param>
    /// <returns></returns>
    public double SeparationSqured(XYZ vec1, XYZ vec2, out XYZ dxdydz)
    {
        dxdydz = vec1 - vec2;

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
    /// Учёт периодических граничных условий.
    /// </summary>
    /// <param name="pos"></param>
    /// <returns></returns>
    public XYZ Periodic(XYZ pos)
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

    // public XYZ Periodic(XYZ pos, XYZ p)
    // {
    //     var newPos = XYZ.Zero;
    //     // Ось X.               
    //     if (pos.X > BoxSize)
    //     {
    //         newPos.X = pos.X - BoxSize;
    //         Flux.X += p.X;
    //     }
    //     else if (pos.X < 0)
    //     {
    //         newPos.X = pos.X + BoxSize;
    //         Flux.X -= p.X;
    //     }
    //     else newPos.X = pos.X;
    //
    //     // Ось Y.
    //     if (pos.Y > BoxSize)
    //     {
    //         newPos.Y = pos.Y - BoxSize;
    //         Flux.Y += p.Y;
    //     }
    //     else if (pos.Y < 0)
    //     {
    //         newPos.Y = pos.Y + BoxSize;
    //         Flux.Y -= p.Y;
    //     }
    //     else newPos.Y = pos.Y;
    //
    //     // Ось Z.
    //     if (pos.Z > BoxSize)
    //     {
    //         newPos.Z = pos.Z - BoxSize;
    //         Flux.Z += p.Z;
    //     }
    //     else if (pos.Z < 0)
    //     {
    //         newPos.Z = pos.Z + BoxSize;
    //         Flux.Z -= p.Z;
    //     }
    //     else newPos.Z = pos.Z;
    //
    //     return newPos;
    // }

    // /// <summary>
    // /// Алгоритм Верле.
    // /// </summary>
    // public void Verlet()
    // {
    //     _ke = 0;
    //     _pe = 0;
    //     _virial = 0;
    //
    //     Atoms.ForEach(atom =>
    //     {
    //         var newPos = atom.Velocity * dt + 0.5 * atom.Acceleration * dt * dt;
    //         atom.Position = Periodic(atom.Position + newPos, atom.Velocity * atom.Weight);
    //         atom.PositionNp += newPos;
    //     });
    //
    //     Atoms.ForEach(atom => atom.Velocity += 0.5 * atom.Acceleration * dt);
    //
    //     // Рассчёт новых ускорений.
    //     // Accel();
    //
    //     Atoms.ForEach(atom =>
    //     {
    //         atom.Velocity += 0.5 * atom.Acceleration * dt;
    //         _ke += 0.5 * atom.Velocity.SquaredMagnitude() * atom.Weight;
    //     });
    //
    //     // Рассчёт АКФ скорости.
    //     if (CurrentStep == 1)
    //         _vtList.Clear();
    //     if (_vtList.Count < CountNumberAcf + CountRepeatAcf * StepRepeatAcf)
    //         _vtList.Add(GetVelocitiesAtoms());
    //
    //     CurrentStep++;
    // }

    // /// <summary>
    // /// Вычисление ускорения, pe, virial.
    // /// </summary>
    // private void Accel()
    // {
    //     Atoms.ForEach(atom => atom.Acceleration = XYZ.Zero);
    //
    //     for (var i = 0; i < CountAtoms - 1; i++)
    //     {
    //         var sumForce = XYZ.Zero;
    //         for (var j = i + 1; j < CountAtoms; j++)
    //         {
    //             var rij = SeparationSqured(Atoms[i].Position, Atoms[j].Position, out var dxdydz);
    //
    //             var force = _potential.PotentialDerivative(new object[] { rij, dxdydz });
    //             sumForce += force;
    //             Atoms[i].Acceleration += force / Atoms[i].Weight;
    //             Atoms[j].Acceleration -= force / Atoms[j].Weight;
    //
    //             _pe += (double)_potential.PotentialEnergy(new object[] { rij });
    //         }
    //
    //         _virial += Atoms[i].Position.X * sumForce.X + Atoms[i].Position.Y * sumForce.Y + Atoms[i].Position.Z * sumForce.Z;
    //     }
    // }
}