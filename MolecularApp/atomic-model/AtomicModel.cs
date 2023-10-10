using System;
using System.Collections.Generic;
using MolecularApp.potentials;

namespace MolecularApp.atomic_model;

public partial class AtomicModel
{
    /// <summary>
    /// Список атомов.
    /// </summary>
    public List<Atom> Atoms { get; }
    
    /// <summary>
    /// Получить количество элементов в словаре которые хранят расстояния между атомами (исключает расстояние между соседями для каждого атома).
    /// </summary>
    public int CountDistanceIJ => DistanceBetweenAtoms.Count;
    
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
    public double Ke => _ke / eV;

    private double _ke;

    /// <summary>
    /// Потенциальная энергия системы (эВ).
    /// </summary>
    public double Pe => _pe / eV;

    private double _pe;

    /// <summary>
    /// Полная энергия системы (эВ).
    /// </summary>
    public double Fe => (_pe + _ke) / eV;

    /// <summary>
    /// Температура системы.
    /// </summary>
    public double T => 2 * _ke / (3 * kB * CountAtoms);

    /// <summary>
    /// Давление системы, рассчитанный через вириал (Па).
    /// </summary>
    public double P1 => (_ke + _virial / CountAtoms) / (3 * V);

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
    /// Параметр решётки сплава SnGe.
    /// </summary>
    public double LatticeGeSn => Atom.GetLattice(AtomType.Ge) * FractionGe + Atom.GetLattice(AtomType.Sn) * FractionSn;
    
    /// <summary>
    /// Параметры потенциала для Ge (германий).
    /// </summary>
    private static TersoffPotential.PotentialParams ParamsGe => new(1769, 419.23, 0.31, 0.28, 9.01e-7, 106430, 15.65, 0.75627, -0.43884, 24.451, 17.047);

    /// <summary>
    /// Параметры потенциала для Sn (олово).
    /// </summary>
    private static TersoffPotential.PotentialParams ParamsSn => new(520.4677, 281.4117, 0.34, 0.30, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 15.5, 12.5649);
    // public static TersoffPotential.PotentialParams ParamsSn => new TersoffPotential.PotentialParams(526.46, 296.83, 0.34, 0.30, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 15.3, 12.56);
    // public static TersoffPotential.PotentialParams ParamsSn => new TersoffPotential.PotentialParams(2848, 658.62, 0.32, 0.28, 6.01e-7, 1.4e+5, 14.5, 0.74, -0.502, 22.5, 16.2);

    /// <summary>
    /// Класс потенциала.
    /// </summary>
    private readonly IPotential _potential;

    // Параметры симуляции.
    /// <summary>
    /// Величина временного шага (c).
    /// </summary>
    public double dt { get; set; }

    /// <summary>
    /// Текущий временной шаг.
    /// </summary>
    public int CurrentStep { get; set; }

    // Константы.
    /// <summary>
    /// 1 эВ в Дж.
    /// </summary>
    private const double eV = 1.602176634e-19;

    /// <summary>
    /// Постоянная Больцмана (Дж/К).
    /// </summary>
    private const double kB = 1.380649e-23;

    /// <summary>
    /// Генератор случайных чисел.
    /// </summary>
    private static Random Rnd => new(Guid.NewGuid().GetHashCode());

    /// <summary>
    /// Список начальных координат атомов.
    /// </summary>
    private List<XYZ> _rt0;

    /// <summary>
    /// Список скоростей атомв в разные моменты времени.
    /// </summary>
    private List<List<XYZ>> _vtList;

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

        // Выбор типа решётки.
        CreateDiamondSystem();

        // Выбор типа потенциала.
        _potential = new TersoffPotential(ParamsGe, ParamsSn);

        // Начальная инициализация параметров.
        InitCalculation();

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
            idx = Rnd.Next(0, CountAtoms);
            var atom = Atoms[idx];
            if (atom.Type != AtomType.Sn)
                atom.Type = AtomType.Sn;
        }
    }
}