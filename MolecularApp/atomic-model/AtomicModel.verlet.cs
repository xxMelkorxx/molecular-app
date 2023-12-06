using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace MolecularApp.atomic_model;

public partial class AtomicModel
{
    /// <summary>
    /// Алгоритм Верле для вычисления координат и скоростей атомов на временном шаге.
    /// </summary>
    public void Verlet()
    {
        _virial = 0;

        // Расчёт новых положений атомов.
        Atoms.ForEach(atom =>
        {
            var newPos = atom.Velocity * dt + 0.5 * atom.Acceleration * dt * dt;
            atom.Position = Periodic(atom.Position + newPos, atom.Velocity * atom.Weight);
            atom.PositionNp += newPos;
        });
        // Поиск соседей атомов после их сдвига.
        SearchAtomsNeighbours();
        // Расчёт новых скоростей и ускорений.
        Atoms.ForEach(atom => atom.Velocity += 0.5 * atom.Acceleration * dt);
        Accels();
        Atoms.ForEach(atom => atom.Velocity += 0.5 * atom.Acceleration * dt);

        // Рассчёт АКФ скорости.
        if (CurrentStep == 1)
            _vtList.Clear();
        if (_vtList.Count < CountNumberAcf + CountRepeatAcf * StepRepeatAcf)
            _vtList.Add(GetVelocitiesAtoms());

        CurrentStep++;
    }

    /// <summary>
    /// Объект-заглушка для синхронизации потоков.
    /// </summary>
    private object _taskLocker = new();

    /// <summary>
    /// Поиск соседей для каждого атома системы;
    /// вычисление <seealso cref="DistanceBetweenAtoms"/>.
    /// </summary>
    public void SearchAtomsNeighbours()
    {
        Atoms.ForEach(a => a.Neighbours.Clear());
        DistanceBetweenAtoms.Clear();

        var searchRadius = _potential.GetRadiusCutoff(FisrtFraction) + 1.2 * (1e-3 * SystemLattice);
        var searchRadiusSquared = searchRadius * searchRadius;

        // Расчёт расстояний между атомами системы в пределах (многопоточность).
        Parallel.For(0, CountAtoms, i =>
        {
            var atomI = Atoms[i];
            for (var j = i + 1; j < CountAtoms; j++)
            {
                var atomJ = Atoms[j];
                var distanceSquared = SeparationSqured(atomI.Position, atomJ.Position);
                if (distanceSquared > searchRadiusSquared)
                    continue;

                var radius = Math.Sqrt(distanceSquared);
                var indexes = PairIndexes.GetIndexes(atomI, atomJ);

                lock (_taskLocker)
                {
                    DistanceBetweenAtoms.Add(indexes, radius);
                    atomI.Neighbours.Add(atomJ);
                    atomJ.Neighbours.Add(atomI);
                }
            }
        });

        // Расчёт расстояний между соседями каждого атома системы в пределах (многопоточность);
        Parallel.For(0, CountAtoms, n =>
        {
            var selAtom = Atoms[n];
            for (var i = 0; i < selAtom.Neighbours.Count; i++)
            for (var j = i + 1; j < selAtom.Neighbours.Count; j++)
            {
                var neighI = selAtom.Neighbours[i];
                var neighJ = selAtom.Neighbours[j];
                var indexes = PairIndexes.GetIndexes(neighI, neighJ);
                var distance = Separation(neighI.Position, neighJ.Position);

                lock (_taskLocker)
                    DistanceBetweenAtoms[indexes] = distance;
            }
        });
    }

    /// <summary>
    /// Подсчёт ускорения системы атомов.
    /// </summary>
    public void Accels()
    {
        Parallel.ForEach(Atoms, atom =>
        {
            var force = ThreePointsDiffAtom(atom);
            atom.Acceleration = force / atom.Weight;
            _virial += atom.Position.X * force.X + atom.Position.Y * force.Y + atom.Position.Z * force.Z;
        });
    }

    /// <summary>
    /// Расчёт производной потенциальной энергии по формуле трехточечного дифференцирования.
    /// </summary>
    /// <param name="atom">Выбранный атом.</param>
    /// <param name="eps"></param>
    /// <returns></returns>
    private XYZ ThreePointsDiffAtom(Atom atom, double eps = 1e-3)
    {
        var delta = eps * SystemLattice;
        double potentialLeft, potentialRight;
        XYZ force;

        potentialRight = PotentialEnergyAtomShifted(atom, Periodic(atom.Position + new XYZ(delta, 0, 0)));
        potentialLeft = PotentialEnergyAtomShifted(atom, Periodic(atom.Position - new XYZ(delta, 0, 0)));
        force.X = potentialLeft - potentialRight;

        potentialRight = PotentialEnergyAtomShifted(atom, Periodic(atom.Position + new XYZ(0, delta, 0)));
        potentialLeft = PotentialEnergyAtomShifted(atom, Periodic(atom.Position - new XYZ(0, delta, 0)));
        force.Y = potentialLeft - potentialRight;

        potentialRight = PotentialEnergyAtomShifted(atom, Periodic(atom.Position + new XYZ(0, 0, delta)));
        potentialLeft = PotentialEnergyAtomShifted(atom, Periodic(atom.Position - new XYZ(0, 0, delta)));
        force.Z = potentialLeft - potentialRight;

        return force / (2 * delta);
    }

    /// <summary>
    /// Подсчёт потенциальной энергии атома и его окружения при сдвиге.
    /// </summary>
    /// <param name="atom">Выбранный атом.</param>
    /// <param name="shiftedPosition">Сдвинутая координата атома.</param>
    /// <returns></returns>
    private double PotentialEnergyAtomShifted(Atom atom, XYZ shiftedPosition)
    {
        var atomsDistances = new Dictionary<PairIndexes, double>();
        atom.Neighbours.ForEach(neigh =>
        {
            var indexes = PairIndexes.GetIndexes(atom, neigh);
            atomsDistances[indexes] = Separation(shiftedPosition, neigh.Position);

            neigh.Neighbours.ForEach(neighK =>
            {
                if (neighK.Index != atom.Index)
                {
                    indexes = PairIndexes.GetIndexes(neigh, neighK);
                    atomsDistances[indexes] = DistanceBetweenAtoms[indexes];
                }
            });
        });

        // Вычисляем расстояние для угловых частей потенциалов.
        RecalcDistanceJK(shiftedPosition, atom, atom, ref atomsDistances);
        foreach (var neigh in atom.Neighbours)
            RecalcDistanceJK(shiftedPosition, atom, neigh, ref atomsDistances);

        return _potential.PotentialDerivative(atom, atomsDistances);
    }

    /// <summary>
    /// Расчёт расстояний между соседями выбранного атома.
    /// </summary>
    /// <param name="shiftedPosition">Сдвинутая координата атома.</param>
    /// <param name="shiftedAtom">Сдвигаемый атом.</param>
    /// <param name="selAtom">Выбранный атом.</param>
    /// <param name="distancesIJ">Словарь, в который записываются расстояния между соседями выбранного атома.</param>
    private void RecalcDistanceJK(XYZ shiftedPosition, Atom shiftedAtom, Atom selAtom, ref Dictionary<PairIndexes, double> distancesIJ)
    {
        for (var i = 0; i < selAtom.Neighbours.Count; i++)
        for (var j = i + 1; j < selAtom.Neighbours.Count; j++)
        {
            var neighI = selAtom.Neighbours[i];
            var neighJ = selAtom.Neighbours[j];
            var indexes = PairIndexes.GetIndexes(neighI, neighJ);

            double distance;
            if (neighI.Index == shiftedAtom.Index)
                distance = Separation(shiftedPosition, neighJ.Position);
            else if (neighJ.Index == shiftedAtom.Index)
                distance = Separation(neighI.Position, shiftedPosition);
            else
                distance = DistanceBetweenAtoms[indexes];

            distancesIJ[indexes] = distance;
        }
    }

    /// <summary>
    /// Вычисление расстояния между частицами с учётом периодических граничных условий. 
    /// </summary>
    /// <param name="vec1"></param>
    /// <param name="vec2"></param>
    /// <returns></returns>
    private double Separation(XYZ vec1, XYZ vec2) => Math.Sqrt(SeparationSqured(vec1, vec2));

    /// <summary>
    /// Вычисление квадрата расстояния между частицами с учётом периодических граничных условий. 
    /// </summary>
    /// <param name="vec1"></param>
    /// <param name="vec2"></param>
    /// <returns></returns>
    public double SeparationSqured(XYZ vec1, XYZ vec2)
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

    public XYZ Periodic(XYZ pos, XYZ p)
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
}