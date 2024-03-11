using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MolecularApp.potentials;

namespace MolecularApp.atomic_model;

public class MonocrystalModel : AtomicModel
{
    private MLJPotential _potential; 
    
    public List<XYZ> rt0;
    
    public AtomType AtomType;
    
    //  Создание атомной модели монокристалла.
    public MonocrystalModel(int size, AtomType atomType)
    {
        Atoms = new List<Atom>();
        Size = size;
        AtomType = atomType;
        CurrentStep = 1;
        Flux = XYZ.Zero;

        // Вычисление параметра решётки системы по закону Вегарда.
        SystemLattice = Atom.GetLattice(atomType);

        _virial = 0;
        _rnd = new Random(Guid.NewGuid().GetHashCode());

        // Инициализация потенциала.
        _potential = new MLJPotential(atomType);

        // Получение начальных координат без учёта ПГУ для первого типа атома и для второго и общий.
        _vtList = new List<List<XYZ>> { GetVelocitiesAtoms() };
    }
    
    public override string GetNameLogFile()
    {
        return $"Results_{DateTime.Now:ddmmyyyy_hhmmss}_{AtomType}_{CountAtoms}";
    }
    
    public override void CreateSystem()
    {
        var idx = 0;
        // Размещение атомов в ГЦК-решётку.
        for (var i = 0; i < Size; i++)
        for (var j = 0; j < Size; j++)
        for (var k = 0; k < Size; k++)
        {
            Atoms.Add(new Atom(++idx, AtomType, new XYZ(i, j, k) * SystemLattice));
            Atoms.Add(new Atom(++idx, AtomType, new XYZ(i + 0.5, j, k + 0.5) * SystemLattice));
            Atoms.Add(new Atom(++idx, AtomType, new XYZ(i, j + 0.5, k + 0.5) * SystemLattice));
            Atoms.Add(new Atom(++idx, AtomType, new XYZ(i + 0.5, j + 0.5, k) * SystemLattice));
        }
    }

    public override void InitCalculation()
    {
        _pe = 0;
        _virial = 0;
        
        Accels();
        
        _ke = Atoms.Sum(atom => atom.Weight * atom.Velocity.SquaredMagnitude() / 2d);
    }

    public override double GetRadiusAtom()
    {
        return _potential.Sigma;
    }

    public override void Verlet()
    {
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
        Accels();

        Atoms.ForEach(atom =>
        {
            atom.Velocity += 0.5 * atom.Acceleration * dt;
            _ke += 0.5 * atom.Velocity.SquaredMagnitude() * atom.Weight;
        });
        
        _ke = Atoms.Sum(atom => atom.Weight * atom.Velocity.SquaredMagnitude() / 2d);

        // Рассчёт АКФ скорости.
        if (CurrentStep == 1)
            _vtList.Clear();
        if (_vtList.Count < CountNumberAcf + CountRepeatAcf * StepRepeatAcf)
            _vtList.Add(GetVelocitiesAtoms());

        CurrentStep++;
    }

    private void Accels()
    {
        Atoms.ForEach(atom => atom.Acceleration = XYZ.Zero);

        for (var i = 0; i < CountAtoms - 1; i++)
        {
            var sumForce = XYZ.Zero;
            for (var j = i + 1; j < CountAtoms; j++)
            {
                var rij = SeparationSqured(Atoms[i].Position, Atoms[j].Position, out var dxdydz);
        
                var force = _potential.PotentialDerivative(rij, dxdydz);
                sumForce += force;
                Atoms[i].Acceleration += force / Atoms[i].Weight;
                Atoms[j].Acceleration -= force / Atoms[i].Weight;
        
                _pe += _potential.PotentialEnergy(rij);
            }
            _virial += Atoms[i].Position.X * sumForce.X + Atoms[i].Position.Y * sumForce.Y + Atoms[i].Position.Z * sumForce.Z;
        }
    }
    
    /// <summary>
    /// Вычисление квадрата расстояния между частицами с учётом периодических граничных условий. 
    /// </summary>
    /// <param name="vec1"></param>
    /// <param name="vec2"></param>
    /// <param name="dxdydz"></param>
    /// <returns></returns>
    private double SeparationSqured(XYZ vec1, XYZ vec2, out XYZ dxdydz)
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
    /// Получение координат атомов без учёта ПГУ на текущем шаге.
    /// </summary>
    public List<XYZ> GetPosNpAtoms() => Atoms.Select(atom => atom.PositionNp).ToList();

    /// <summary>
    /// Вычисление среднего квадрата смещения на текущем шаге. 
    /// </summary>
    public double GetMsd() => rt0.Zip(GetPosNpAtoms(), (vec1, vec2) => (vec2 - vec1).SquaredMagnitude()).Sum() / CountAtoms;
}