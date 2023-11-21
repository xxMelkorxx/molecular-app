using System;
using System.Collections.Generic;
using System.Drawing;
using ScottPlot;
using ScottPlot.Control;

namespace MolecularApp;

public partial class MainWindow
{
    private class FindPrimesInput
    {
        public List<object> Args;

        public FindPrimesInput(IEnumerable<object> args)
        {
            Args = new List<object>();
            foreach (var arg in args)
                Args.Add(arg);
        }
    }

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
               $"Сплав SnGe: {_atomic.SecondFraction * 100} %:{_atomic.FisrtFraction * 100} %\n" +
               $"Размер структуры (Nx/Ny/Nz) - {_atomic.Size}/{_atomic.Size}/{_atomic.Size}\n" +
               $"Размер структуры (Lx/Ly/Lz) - {(_atomic.BoxSize * 1e9).ToString("F3")}/{(_atomic.BoxSize * 1e9).ToString("F3")}/{(_atomic.BoxSize * 1e9).ToString("F3")} нм\n" +
               $"Объём - {(_atomic.V * 1e27).ToString("F5")} нм³\n" +
               $"Число атомов - {_atomic.CountAtoms}\n" +
               $"Параметр решётки - {(_atomic.SystemLattice * 1e9).ToString("F3")} нм\n" +
               $"Кинетическая энергия - {_atomic.Ke.ToString("F5")} эВ\n" +
               $"Потенциальная энергия - {_atomic.Pe.ToString("F5")} эВ\n" +
               $"Полная энергия - {_atomic.Fe.ToString("F5")} эВ\n";
        // $"Температура - {_atomic.T.ToString("F1")} К\n" +
        // $"Давление - {Math.Round(_atomic.P1)} Па\n";
    }

    /// <summary>
    /// Заголовок таблицы.
    /// </summary>
    /// <returns></returns>
    private static string TableHeader()
    {
        return $"{"Шаг".PadLeft(6)} |" +
               $"{"Кин.энергия(эВ)".PadLeft(16)} |" +
               $"{"Пот.энергия(эВ)".PadLeft(16)} |" +
               $"{"Полн.энергия(эВ)".PadLeft(17)} |\n";
               // + $"{"Температура(К)".PadLeft(15)} |"
               // + $"{"Давление 1(Па)".PadLeft(15)} |"
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
        return $"{i.ToString().PadLeft(6)} |" +
               $"{_atomic.Ke.ToString("F5").PadLeft(16)} |" +
               $"{_atomic.Pe.ToString("F5").PadLeft(16)} |" +
               $"{_atomic.Fe.ToString("F5").PadLeft(17)} |\n";
               // + $"{_atomic.T.ToString("F1").PadLeft(15)} |"
               // + $"{Math.Round(_atomic.P1).ToString("F0").PadLeft(15)} |\n";
               // + $"{(_atomic.P2 / nsnap).ToString("F1").PadLeft(15)} |\n";
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