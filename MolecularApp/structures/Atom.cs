using System;
using System.Collections.Generic;

namespace MolecularApp;

public enum AtomType
{
    Sn,
    Ge
}

public class Atom
{
    /// <summary>
    /// Идентификатор.
    /// </summary>
    public int Index { get; set; }

    /// <summary>
    /// Вектор координаты.
    /// </summary>
    public XYZ Position { get; set; }

    /// <summary>
    /// Координаты атома без учёта периодичности границ.
    /// </summary>
    public XYZ PositionNp { get; set; }

    /// <summary>
    /// Вектор скорости.
    /// </summary>
    public XYZ Velocity { get; set; }

    /// <summary>
    /// Вектор ускорения.
    /// </summary>
    public XYZ Acceleration { get; set; }

    /// <summary>
    /// Масса атома (кг).
    /// </summary>
    public double Weight => Type switch
    {
        AtomType.Sn => 118.71 * 1.66054e-27,
        AtomType.Ge => 72.63 * 1.66054e-27,
        _ => throw new ArgumentNullException()
    };

    /// <summary>
    /// Тип атома.
    /// </summary>
    public AtomType Type { get; set; }

    /// <summary>
    /// Список соседей атома.
    /// </summary>   
    public List<Atom> Neighbours;

    /// <summary>
    /// Создание атома.
    /// </summary>
    /// <param name="index"></param>
    /// <param name="atomType"></param>
    /// <param name="pos"></param>
    public Atom(int index, AtomType atomType, XYZ pos)
    {
        Index = index;
        Type = atomType;
        Position = pos;
        PositionNp = pos;
        Velocity = XYZ.Zero;
        Acceleration = XYZ.Zero;
        Neighbours = new List<Atom>();
    }

    /// <summary>
    /// Получение массы атома.
    /// </summary>
    /// <param name="type"></param>
    /// <returns></returns>
    /// <exception cref="ArgumentNullException"></exception>
    public static double GetWeightAtom(AtomType type) => type switch
    {
        AtomType.Sn => 118.71 * 1.66054e-27,
        AtomType.Ge => 72.63 * 1.66054e-27
    };

    /// <summary>
    /// Параметр решётки (м)
    /// </summary>
    public static double GetLattice(AtomType type) => type switch
    {
        AtomType.Sn => 0.54307e-9,
        AtomType.Ge => 0.566e-9
    };
}