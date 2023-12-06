using System;
using System.Collections.Generic;
using System.Linq;
using MolecularApp.potentials;

namespace MolecularApp.atomic_model;

public partial class AtomicModel
{
    /// <summary>
    /// Список атомов.
    /// </summary>
    public List<Atom> Atoms { get; }

    /// <summary>
    /// Расстояния между атомами с учётом параметра обрезания выбранного потенциала.
    /// </summary>
    public Dictionary<PairIndexes, double> DistanceBetweenAtoms { get; }

    /// <summary>
    /// Размер расчётной ячейки.
    /// </summary>
    public int Size { get; }

    /// <summary>
    /// Размер расчётной ячейки (м).
    /// </summary>
    public double BoxSize => Size * SystemLattice;

    /// <summary>
    /// Число атомов.
    /// </summary>
    public int CountAtoms => Atoms.Count;

    /// <summary>
    /// Кинетическая энергия системы (эВ).
    /// </summary>
    public double Ke => Atoms.Sum(atom => atom.Weight * atom.Velocity.SquaredMagnitude() / 2d) / eV;

    /// <summary>
    /// Потенциальная энергия системы (эВ).
    /// </summary>
    public double Pe => Atoms.Sum(atom => _potential.PotentialEnergy(atom, DistanceBetweenAtoms)) / eV;

    /// <summary>
    /// Полная энергия системы (эВ).
    /// </summary>
    public double Fe => (Ke + Pe);

    /// <summary>
    /// Температура системы.
    /// </summary>
    public double T => 2 * Ke * eV / (3 * kB * CountAtoms);

    /// <summary>
    /// Давление системы, рассчитанный через вириал (Па).
    /// </summary>
    public double P1 => (Ke * eV + _virial / CountAtoms) / (3 * V);

    private double _virial;
    
    /// <summary>
    /// Давление системы (Па).
    /// </summary>
    public double P2 => (Flux.X + Flux.Y + Flux.Z) / (6d * BoxSize * BoxSize * dt);

    public XYZ Flux;

    /// <summary>
    /// Объём системы (м³).
    /// </summary>
    public double V => BoxSize * BoxSize * BoxSize;

    /// <summary>
    /// Доля первого элемента в сплаве.
    /// </summary>
    public double FisrtFraction { get; set; }

    /// <summary>
    /// Доля второго элемента в сплаве.
    /// </summary>
    public double SecondFraction { get; set; }

    /// <summary>
    /// Параметр решётки сплава (м).
    /// </summary>
    public double SystemLattice { get; set; }

    // Параметры симуляции.
    /// <summary>
    /// Величина временного шага (c).
    /// </summary>
    public double dt { get; set; }

    /// <summary>
    /// Текущий временной шаг.
    /// </summary>
    public int CurrentStep { get; set; }

    /// <summary>
    /// Число отсчётов.
    /// </summary>
    public int CountNumberAcf;

    /// <summary>
    /// Число повторений.
    /// </summary>
    public int CountRepeatAcf;

    /// <summary>
    /// Шаг повторений подсчёта.
    /// </summary>
    public int StepRepeatAcf;

    // Константы.
    private const double eV = 1.602176634e-19; // 1 эВ в Дж.
    private const double kB = 1.380649e-23; // Постоянная Больцмана (Дж/К).

    /// <summary>
    /// Класс потенциала.
    /// </summary>
    private readonly IPotential _potential;

    private readonly Random _rnd;
    private List<XYZ> _rt0;
    private List<List<XYZ>> _vtList;

    /// <summary>
    /// Создание атомной модели.
    /// </summary>
    /// <param name="size"></param>
    /// <param name="firstTypeAtom"></param>
    /// <param name="fisrtFraction"></param>
    /// <param name="secondTypeAtom"></param>
    /// <param name="secondFraction"></param>
    public AtomicModel(int size, AtomType firstTypeAtom, double fisrtFraction, AtomType secondTypeAtom, double secondFraction)
    {
        Atoms = new List<Atom>();
        DistanceBetweenAtoms = new Dictionary<PairIndexes, double>();
        Size = size;
        FisrtFraction = fisrtFraction;
        SecondFraction = secondFraction;
        CurrentStep = 1;
        Flux = XYZ.Zero;

        // Вычисление параметра решётки системы по закону Вегарда.
        SystemLattice = Atom.GetLattice(firstTypeAtom) * FisrtFraction + Atom.GetLattice(secondTypeAtom) * SecondFraction;

        _virial = 0;
        _rnd = new Random(Guid.NewGuid().GetHashCode());

        // Инициализация потенциала.
        _potential = new TersoffPotential(firstTypeAtom, secondTypeAtom);
        // Создание расчётной ячейки.
        CreateDiamondSystem(firstTypeAtom, secondTypeAtom);
        // Начальный расчёт характеристик.
        InitCalculation();

        _rt0 = GetPosNpAtoms();
        _vtList = new List<List<XYZ>> { GetVelocitiesAtoms() };
    }

    /// <summary>
    /// Начальное размещение атомов в АЦК-решётку.
    /// </summary>
    private void CreateDiamondSystem(AtomType firstTypeAtom, AtomType secondTypeAtom)
    {
        var idx = 0;
        for (var i = 0; i < Size; i++)
        for (var j = 0; j < Size; j++)
        for (var k = 0; k < Size; k++)
        {
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i, j, k) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i + 0.5, j, k + 0.5) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i, j + 0.5, k + 0.5) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i + 0.5, j + 0.5, k) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i + 0.25, j + 0.25, k + 0.25) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i + 0.25, j + 0.75, k + 0.75) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i + 0.75, j + 0.25, k + 0.75) * SystemLattice));
            Atoms.Add(new Atom(++idx, firstTypeAtom, new XYZ(i + 0.75, j + 0.75, k + 0.25) * SystemLattice));
        }

        // Заполнение системы атомами олова.
        var countSwapAtoms = (int)(CountAtoms * SecondFraction);
        for (var i = 0; i < countSwapAtoms; i++)
        {
            idx = _rnd.Next(0, CountAtoms);
            if (Atoms[idx].Type != secondTypeAtom)
                Atoms[idx].Type = secondTypeAtom;
            else i--;
        }
    }

    /// <summary>
    /// Вычисление начальных параметров системы (Acceleration, Pe, Ke, Press).
    /// </summary>
    public void InitCalculation()
    {
        // Рассчёт соседей для каждого атома.
        SearchAtomsNeighbours();
        // Начальный подсчёт ускорений атомов.
        Accels();
    }
}