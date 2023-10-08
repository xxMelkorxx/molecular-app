namespace MolecularApp;

public struct PairIndexes
{
    public int Index1 { get; set; }
    public int Index2 { get; set; }

    public PairIndexes(int idx1, int idx2)
    {
        Index1 = idx1;
        Index2 = idx2;
    }

    /// <summary>
    /// Получить структуру упорядоченных индексов двух атомов.
    /// </summary>
    /// <param name="atomI">Первый атом.</param>
    /// <param name="atomJ">Второй атом.</param>
    /// <returns></returns>
    public static PairIndexes GetIndexes(Atom atomI, Atom atomJ) => atomI.Index < atomJ.Index ? new(atomI.Index, atomJ.Index) : new(atomJ.Index, atomI.Index);
}