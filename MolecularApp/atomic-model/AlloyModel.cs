using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MolecularApp.potentials;

namespace MolecularApp.atomic_model;

public class AlloyModel : AtomicModel
{
    // Типы атомов.
    public AtomType FirstAtomType, SecondAtomType;
    
    // Доля первого элемента в сплаве.
    public double FisrtFraction { get; }
    // Доля второго элемента в сплаве.
    public double SecondFraction { get; }
    
    // Расстояния между атомами с учётом параметра обрезания выбранного потенциала.
    private Dictionary<PairIndexes, double> DistanceBetweenAtoms { get; }

    //  Создание атомной модели сплава.
    public AlloyModel(int size, AtomType firstTypeAtom, double fisrtFraction, AtomType secondTypeAtom, double secondFraction)
    {
        Atoms = new List<Atom>();
        Size = size;
        FirstAtomType = firstTypeAtom;
        SecondAtomType = secondTypeAtom;
        FisrtFraction = fisrtFraction;
        SecondFraction = secondFraction;
        CurrentStep = 1;
        Flux = XYZ.Zero;

        // Вычисление параметра решётки системы по закону Вегарда.
        SystemLattice = Atom.GetLattice(firstTypeAtom) * fisrtFraction + Atom.GetLattice(secondTypeAtom) * secondFraction;

        _virial = 0;
        _rnd = new Random(Guid.NewGuid().GetHashCode());

        // Инициализация потенциала.
        _potential = new TersoffPotential(firstTypeAtom, secondTypeAtom);
        DistanceBetweenAtoms = ((TersoffPotential)_potential).AtomsDistances;

        // Получение начальных координат без учёта ПГУ для первого типа атома и для второго и общий.
        _vtList = new List<List<XYZ>> { GetVelocitiesAtoms() };
    }
    
    public override void CreateSystem()
    {
        var idx = 0;
        // Размещение атомов в АЦК-решётку.
        for (var i = 0; i < Size; i++)
        for (var j = 0; j < Size; j++)
        for (var k = 0; k < Size; k++)
        {
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i, j, k) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i + 0.5, j, k + 0.5) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i, j + 0.5, k + 0.5) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i + 0.5, j + 0.5, k) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i + 0.25, j + 0.25, k + 0.25) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i + 0.25, j + 0.75, k + 0.75) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i + 0.75, j + 0.25, k + 0.75) * SystemLattice));
            Atoms.Add(new Atom(++idx, FirstAtomType, new XYZ(i + 0.75, j + 0.75, k + 0.25) * SystemLattice));
        }

        // Заполнение системы атомами олова.
        var countSwapAtoms = (int)(CountAtoms * SecondFraction);
        for (var i = 0; i < countSwapAtoms; i++)
        {
            idx = _rnd.Next(0, CountAtoms);
            if (Atoms[idx].Type != SecondAtomType)
                Atoms[idx].Type = SecondAtomType;
            else i--;
        }
    }
    
    public override void InitCalculation()
    {
        // Рассчёт соседей для каждого атома.
        SearchAtomsNeighbours();
        // Начальный подсчёт ускорений атомов.
        Accels();
    }

    public override double GetRadiusAtom() => ((TersoffPotential)_potential).GetRadiusCutoff(FisrtFraction);
    
    public override void Verlet()
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

    protected override void Accels()
    {
        Parallel.ForEach(Atoms, atom =>
        {
            var force = ThreePointsDiffAtom(atom);
            atom.Acceleration = force / atom.Weight;
            _virial += atom.Position.X * force.X + atom.Position.Y * force.Y + atom.Position.Z * force.Z;
        });
    }

    // Объект-заглушка для синхронизации потоков.
    private object _taskLocker = new();

    // Поиск соседей для каждого атома системы;
    public void SearchAtomsNeighbours()
    {
        Atoms.ForEach(a => a.Neighbours.Clear());
        DistanceBetweenAtoms.Clear();

        var searchRadius = GetRadiusAtom() + 1.2 * SystemLattice * 1e-3;
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
    
    // Расчёт производной потенциальной энергии по формуле трехточечного дифференцирования.
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

    // Подсчёт потенциальной энергии атома и его окружения при сдвиге.
    private double PotentialEnergyAtomShifted(Atom atom, XYZ shiftedPosition)
    {
        atom.Neighbours.ForEach(neigh =>
        {
            var indexes = PairIndexes.GetIndexes(atom, neigh);
            DistanceBetweenAtoms[indexes] = Separation(shiftedPosition, neigh.Position);
        });

        // Вычисляем расстояние для угловых частей потенциалов.
        RecalcDistanceJK(shiftedPosition, atom, atom);
        foreach (var neigh in atom.Neighbours)
            RecalcDistanceJK(shiftedPosition, atom, neigh);

        return _potential.PotentialDerivative(atom);
    }

    /// <summary>
    /// Расчёт расстояний между соседями выбранного атома.
    /// </summary>
    /// <param name="shiftedPosition">Сдвинутая координата атома.</param>
    /// <param name="shiftedAtom">Сдвигаемый атом.</param>
    /// <param name="selAtom">Выбранный атом.</param>
    private void RecalcDistanceJK(XYZ shiftedPosition, Atom shiftedAtom, Atom selAtom)
    {
        for (var i = 0; i < selAtom.Neighbours.Count; i++)
        for (var j = i + 1; j < selAtom.Neighbours.Count; j++)
        {
            var neighI = selAtom.Neighbours[i];
            var neighJ = selAtom.Neighbours[j];
            var indexes = PairIndexes.GetIndexes(neighI, neighJ);
            
            if (neighI.Index == shiftedAtom.Index)
                DistanceBetweenAtoms[indexes] = Separation(shiftedPosition, neighJ.Position);
            else if (neighJ.Index == shiftedAtom.Index)
                DistanceBetweenAtoms[indexes] = Separation(neighI.Position, shiftedPosition);
        }
    }
    
    // Получение координат атомов без учёта ПГУ на текущем шаге.
    public List<XYZ> GetPosNpAtoms(int flag = 0) => flag switch
    {
        1 => Atoms.Where(atom => atom.Type == FirstAtomType).Select(atom => atom.PositionNp).ToList(),
        2 => Atoms.Where(atom => atom.Type == SecondAtomType).Select(atom => atom.PositionNp).ToList(),
        _ => GetPosNpAtoms()
    };

    // Вычисление среднего квадрата смещения на текущем шаге.
    public double GetMsd(int flag = 0) => flag switch
    {
        1 => rt01.Zip(GetPosNpAtoms(flag), (vec1, vec2) => (vec2 - vec1).SquaredMagnitude()).Sum() / Atoms.Count(atom => atom.Type == FirstAtomType),
        2 => rt02.Zip(GetPosNpAtoms(flag), (vec1, vec2) => (vec2 - vec1).SquaredMagnitude()).Sum() / Atoms.Count(atom => atom.Type == SecondAtomType),
        _ => GetMsd()
    };
}