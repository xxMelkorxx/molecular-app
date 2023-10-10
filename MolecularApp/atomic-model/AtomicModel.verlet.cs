using System;

namespace MolecularApp.atomic_model;

public partial class AtomicModel
{
    /// <summary>
    /// Вычисление начальных параметров системы (Acceleration, Pe, Ke, Press).
    /// </summary>
    public void InitCalculation()
    {
        _ke = 0;
        _pe = 0;
        _virial = 0;
        
        // Рассчёт радиуса поиска соседей с учётом величины сдвига при расчёте ускорений (методом трёхточечного дифференцирования).
        var searchRadius = _potential.GetRadiusCutoff(FractionGe) + 1.2 * (1e-3 * LatticeGeSn);
        
        // Рассчёт соседей для каждого атома.
        SearchAtomsNeighbours(searchRadius);
        
        Accel();

        Atoms.ForEach(atom => _ke += 0.5 * atom.Velocity.SquaredMagnitude() * atom.Weight);
    }

    /// <summary>
    /// Алгоритм Верле.
    /// </summary>
    public void Verlet()
    {
        _ke = 0;
        _pe = 0;
        _virial = 0;

        Atoms.ForEach(atom =>
        {
            var newPos = atom.Velocity * dt + 0.5 * atom.Acceleration * dt * dt;
            atom.Position = Periodic(atom.Position + newPos, atom.Velocity * atom.Weight);
            atom.PositionNp += newPos;
        });

        Atoms.ForEach(atom => atom.Velocity += 0.5 * atom.Acceleration * dt);

        // Рассчёт новых ускорений.
        Accel();

        Atoms.ForEach(atom =>
        {
            atom.Velocity += 0.5 * atom.Acceleration * dt;
            _ke += 0.5 * atom.Velocity.SquaredMagnitude() * atom.Weight;
        });

        // Рассчёт АКФ скорости.
        if (CurrentStep == 1)
            _vtList.Clear();
        if (_vtList.Count < CountNumberAcf + CountRepeatAcf * StepRepeatAcf)
            _vtList.Add(GetVelocitiesAtoms());

        CurrentStep++;
    }

    /// <summary>
    /// Вычисление ускорения, pe, virial.
    /// </summary>
    private void Accel()
    {
        Atoms.ForEach(atom => atom.Acceleration = XYZ.Zero);

        for (var i = 0; i < CountAtoms - 1; i++)
        {
            var sumForce = XYZ.Zero;
            for (var j = i + 1; j < CountAtoms; j++)
            {
                var rij = SeparationSqured(Atoms[i].Position, Atoms[j].Position, out var dxdydz);

                var force = _potential.PotentialDerivative(new object[] { rij, dxdydz });
                sumForce += force;
                Atoms[i].Acceleration += force / Atoms[i].Weight;
                Atoms[j].Acceleration -= force / Atoms[j].Weight;

                _pe += (double)_potential.PotentialEnergy(new object[] { rij });
            }
            _virial += Atoms[i].Position.X * sumForce.X + Atoms[i].Position.Y * sumForce.Y + Atoms[i].Position.Z * sumForce.Z;
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
    /// <param name="p"></param>
    /// <returns></returns>
    public XYZ Periodic(XYZ pos, XYZ p)
    {
        var newPos = XYZ.Zero;
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