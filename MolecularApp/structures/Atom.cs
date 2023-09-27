namespace MolecularApp;

public enum AtomType
{
    Sn,
    Si,
    Ge
}

public class Atom
{
    public Vector Position;
    public Vector Velocity;
    public Vector Acceleration;
    public AtomType Type;

    /// <summary>
    /// Индекс атома.
    /// </summary>
    public ushort Index { get; private set; }
    
    public Atom(ushort idx, AtomType atomType, Vector pos)
    {
        this.Type = atomType;
        this.Position = pos;
        this.Velocity = Vector.Zero;
        this.Acceleration = Vector.Zero;
        this.Index = idx;
        this.Type = atomType;
    }
}