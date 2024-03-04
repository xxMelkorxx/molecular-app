using System;
using System.Collections.Generic;
using System.Linq;
using MolecularApp.potentials;

namespace MolecularApp.atomic_model;

public class MonocrystalModel : AtomicModel
{
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
        _potential = new LJPotential();

        // Получение начальных координат без учёта ПГУ для первого типа атома и для второго и общий.
        _vtList = new List<List<XYZ>> { GetVelocitiesAtoms() };
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
        
    }

    public override double GetRadiusAtom()
    {
        return 0;
    }

    public override void Verlet()
    {
        
    }

    protected override void Accels()
    {
        
    }
}