using System;
using System.Drawing;
using ScottPlot;
using ScottPlot.Control;

namespace MolecularApp;

public partial class MainWindow
{
    private static void SetUpChart(IPlotControl chart, string title, string labelX, string labelY)
    {
        chart.Plot.Title(title);
        chart.Plot.XLabel(labelX);
        chart.Plot.YLabel(labelY);
        chart.Plot.XAxis.MajorGrid(enable: true, color: Color.FromArgb(50, Color.Black));
        chart.Plot.YAxis.MajorGrid(enable: true, color: Color.FromArgb(50, Color.Black));
        chart.Plot.XAxis.MinorGrid(enable: true, color: Color.FromArgb(30, Color.Black), lineStyle: LineStyle.Dot);
        chart.Plot.YAxis.MinorGrid(enable: true, color: Color.FromArgb(30, Color.Black), lineStyle: LineStyle.Dot);
        chart.Plot.Margins(x: 0.0, y: 0.6);
        chart.Plot.SetAxisLimits(xMin: 0, yMin: 0);
        chart.Refresh();
    }

    /// <summary>
    /// Начальная информация о системе.
    /// </summary>
    /// <returns></returns>
    private string InitInfoSystem()
    {
        return "Структура создана...\n" +
               $"Сплав SnGe: {_atomic.FisrtFraction * 100}%:{_atomic.SecondFraction * 100}%\n" +
               $"Размер структуры (Nx/Ny/Nz) - {_atomic.Size}/{_atomic.Size}/{_atomic.Size}\n" +
               $"Размер структуры (Lx/Ly/Lz) - {_atomic.BoxSize * 1e9:F3}/{_atomic.BoxSize * 1e9:F3}/{_atomic.BoxSize * 1e9:F3} нм\n" +
               $"Объём - {_atomic.GetVolume * 1e27:F5} нм³\n" +
               $"Число атомов - {_atomic.CountAtoms}\n" +
               $"Параметр решётки - {_atomic.SystemLattice * 1e9:F3} нм\n" +
               $"Кинетическая энергия - {_atomic.Ke:F5} эВ\n" +
               $"Потенциальная энергия - {_atomic.Pe:F5} эВ\n" +
               $"Полная энергия - {_atomic.Fe:F5} эВ\n";
    }

    /// <summary>
    /// Заголовок таблицы.
    /// </summary>
    /// <returns></returns>
    private static string TableHeader()
    {
        return $"{"Шаг",6} |" +
               $"{"Кин.энергия(эВ)",16} |" +
               $"{"Пот.энергия(эВ)",16} |" +
               $"{"Полн.энергия(эВ)",17} |" +
               $"{"Температура(К)",15} |" +
               $"{"Давление (МПа)",15} |" +
               $"{"Объём(нм³)",11} |\n";
        // + $"{"Давление 2(Па)".PadLeft(15)} |\n";
    }

    /// <summary>
    /// Вывод данных в таблицу.
    /// </summary>
    /// <param name="i"></param>
    /// <param name="nsnap"></param>
    /// <returns></returns>
    private string TableData(int i, int nsnap)
    {
        return $"{i,6} |" +
               $"{_atomic.Ke,16:F5} |" +
               $"{_atomic.Pe,16:F5} |" +
               $"{_atomic.Fe,17:F5} |" +
               $"{_atomic.T,15:F1} |" +
               $"{_atomic.P1 / 1e6,15:F1} |" +
               $"{_atomic.GetVolume * 1e27,11:F5} |\n";
    }

    /// <summary>
    /// Звуковой сигнал.
    /// </summary>
    /// <param name="freq">Частота сигнала.</param>
    /// <param name="duration">Длительность сигнала (мкс).</param>
    /// <param name="count">Повторений.</param>
    private static void AlarmBeep(int freq, int duration, int count)
    {
        for (var i = 0; i < count; i++)
            Console.Beep(freq, duration);
    }
}