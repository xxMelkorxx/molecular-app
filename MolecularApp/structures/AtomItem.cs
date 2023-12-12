namespace MolecularApp;

public struct AtomItem
{
    public AtomType Type;
    public XYZ Position;

    public AtomItem(AtomType type, XYZ position)
    {
        Type = type;
        Position = position;
    }
}