using System;
using System.Collections.Generic;
using System.Linq;
using MolecularApp.potentials;

namespace MolecularApp.atomic_model;

public abstract partial class AtomicModel
{
    // Константы.
    protected const double eV = 1.602176634e-19; // 1 эВ в Дж.
    protected const double kB = 1.380649e-23; // Постоянная Больцмана (Дж/К).
    protected Random _rnd;
    protected List<List<XYZ>> _vtList;
    
    // Список атомов.
    public List<Atom> Atoms { get; protected set; }
    
    // Число атомов.
    public int CountAtoms => Atoms.Count;
    
    // Размер расчётной ячейки.
    public int Size { get; set; }
        
    // Параметр решётки сплава (м).
    public double SystemLattice { get; set; }

    // Размер расчётной ячейки (м).
    public double BoxSize => Size * SystemLattice;
    
    // Кинетическая энергия системы (эВ).
    public double Ke => _ke / eV;
    protected double _ke;

    // Потенциальная энергия системы (эВ).
    public double Pe => _pe / eV;
    protected double _pe;
    
    // Полная энергия системы (эВ).
    public double Fe => Ke + Pe;
    
    // Температура системы.
    public double T => 2 * Ke * eV / (3 * kB * CountAtoms);
    
    // Давление системы, рассчитанный через вириал (Па).
    public double P1 => (Ke * eV + _virial / CountAtoms) / (3 * GetVolume);
    protected double _virial;
    
    // Давление системы (Па).
    public double P2 => (Flux.X + Flux.Y + Flux.Z) / (6d * BoxSize * BoxSize * dt);

    public XYZ Flux;
    
    // Объём системы (м³).
    public double GetVolume => BoxSize * BoxSize * BoxSize;
    
    // Величина временного шага (c).
    public double dt;
    // Текущий временной шаг.
    public int CurrentStep;
    // Число отсчётов.
    public int CountNumberAcf;
    // Число повторений.
    public int CountRepeatAcf;
    // Шаг повторений подсчёта.
    public int StepRepeatAcf;

    /// <summary>
    /// Получение наименование для лог-файла.
    /// </summary>
    public abstract string GetNameLogFile();
    
    /// <summary>
    /// Создание атомой системы.
    /// </summary>
    public abstract void CreateSystem();

    /// <summary>
    /// Инициализация начальный параметров.
    /// </summary>
    public abstract void InitCalculation();

    /// <summary>
    /// Получение радиуса атома.
    /// </summary>
    /// <returns></returns>
    public abstract double GetRadiusAtom();

    /// <summary>
    /// Алгоритм Верле для вычисления координат и скоростей атомов на временном шаге.
    /// </summary>
    public abstract void Verlet();
}