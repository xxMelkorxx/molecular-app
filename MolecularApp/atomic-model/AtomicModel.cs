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
    public double BoxSize => Size * LatticeGeSn;

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
    public double T => 2 * Ke / (3 * kB * CountAtoms);

    /// <summary>
    /// Давление системы, рассчитанный через вириал (Па).
    /// </summary>
    public double P1 => (Ke + _virial / CountAtoms) / (3 * V);

    private double _virial;

    /// <summary>
    /// Давление системы (Па).
    /// </summary>
    public double P2 => (Flux.X + Flux.Y + Flux.Z) / (6 * BoxSize * BoxSize * dt);

    public XYZ Flux;

    /// <summary>
    /// Объём системы (м³).
    /// </summary>
    public double V => BoxSize * BoxSize * BoxSize;

    /// <summary>
    /// Доля Ge.
    /// </summary>
    public double FractionGe { get; }

    /// <summary>
    /// Доля Sn.
    /// </summary>
    public double FractionSn => 1 - FractionGe;

    /// <summary>
    /// Параметр решётки сплава SnGe (м).
    /// </summary>
    public double LatticeGeSn => Atom.GetLattice(AtomType.Ge) * FractionGe + Atom.GetLattice(AtomType.Sn) * FractionSn;

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
    /// Параметры потенциала для Ge (германий).
    /// </summary>
    // private static TersoffParams ParamsGe => new(1769, 419.23, 0.31, 0.28, 9.01e-7, 106430, 15.65, 0.75627, -0.43884, 24.451, 17.047);
    private static TersoffParams ParamsGe => new(1769, 419.23, 0.24451, 0.17047, 9.01e-7, 0.75627, 1.0643e+5, 15.65, -0.43884, 2.8, 3.1);

    /// <summary>
    /// Параметры потенциала для Sn (олово).
    /// </summary>
    // private static TersoffParams ParamsSn => new(520.4677, 281.4117, 0.34, 0.30, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 15.5, 12.5649);
    private static TersoffParams ParamsSn => new(520.4677, 281.4117, 0.155, 0.125649, 6.01e-7, 0.74, 1.4e+5, 14.5, -0.502, 3.0, 3.4);
    // public static TersoffPotential.PotentialParams ParamsSn => new TersoffPotential.PotentialParams(526.46, 296.83, 0.34, 0.30, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 15.3, 12.56);
    // public static TersoffPotential.PotentialParams ParamsSn => new TersoffPotential.PotentialParams(2848, 658.62, 0.32, 0.28, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 22.5, 16.2);

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
    /// <param name="fraction"></param>
    public AtomicModel(int size, double fraction)
    {
        Atoms = new List<Atom>();
        DistanceBetweenAtoms = new Dictionary<PairIndexes, double>();
        Size = size;
        FractionGe = fraction;
        CurrentStep = 1;

        _rnd = new Random(Guid.NewGuid().GetHashCode());
        _potential = new TersoffPotential(ParamsGe, ParamsSn); // Инициализация потенциала.
        CreateDiamondSystem(); // Создание расчётной ячейки.
        InitCalculation(); // Начальный расчёт характеристик.

        _rt0 = GetPosNpAtoms();
        _vtList = new List<List<XYZ>> { GetVelocitiesAtoms() };
    }

    /// <summary>
    /// Начальное размещение атомов в АЦК-решётку.
    /// </summary>
    private void CreateDiamondSystem()
    {
        var idx = 0;
        for (var i = 0; i < Size; i++)
        for (var j = 0; j < Size; j++)
        for (var k = 0; k < Size; k++)
        {
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i, j, k) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i + 0.5, j, k + 0.5) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i, j + 0.5, k + 0.5) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i + 0.5, j + 0.5, k) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i + 0.25, j + 0.25, k + 0.25) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i + 0.25, j + 0.75, k + 0.75) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i + 0.75, j + 0.25, k + 0.75) * LatticeGeSn));
            Atoms.Add(new Atom(++idx, AtomType.Ge, new XYZ(i + 0.75, j + 0.75, k + 0.25) * LatticeGeSn));
        }

        // Заполнение системы атомами олова.
        var countSwapAtoms = (int)(CountAtoms * FractionSn);
        for (var i = 0; i < countSwapAtoms; i++)
        {
            idx = _rnd.Next(0, CountAtoms);
            if (Atoms[idx].Type != AtomType.Sn)
                Atoms[idx].Type = AtomType.Sn;
            else i--;
        }
    }
}