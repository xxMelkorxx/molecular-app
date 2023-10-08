using ScottPlot;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Drawing;
using System.Windows;
using System.Windows.Threading;
using MolecularApp.atomic_model;
using MolecularApp.potentials;
using MolecularApp.scene_manager;

namespace MolecularApp;

public partial class MainWindow
{
    private AtomicModel _atomic;
    private SceneManager _scene;
    private readonly BackgroundWorker _bgWorkerCreateModel, _bgWorkerCalculation;
    private readonly System.Windows.Forms.Timer _timer;
    private List<List<XYZ>> _positionsAtomsList;

    private List<PointD> _msdPoints;
    private List<double> _keValues, _peValues, _feValues;
    private double _averT, _averP;
    private bool _isDisplacement, _isSnapshot, _isNormSpeeds, _isNewSystem;
    private double _yMaxRb;
    private int _initStep;

    public MainWindow()
    {
        InitializeComponent();

        _bgWorkerCreateModel = (BackgroundWorker)FindResource("BackgroundWorkerCreateModel");
        _bgWorkerCalculation = (BackgroundWorker)FindResource("BackgroundWorkerCalculation");

        _timer = new() { Interval = 30 };
        _timer.Tick += OnTickTimer;
    }

    private void OnLoadedMainWindow(object sender, RoutedEventArgs e)
    {
        // Настройка графиков.
        SetUpChart(ChartEnergy, "Графики энергий системы", "t, пс", "E, эВ");
        SetUpChart(ChartRadDist, "График радиального распределения системы", "r, нм", "g(r)");
        SetUpChart(ChartMsd, "График  среднего квадрата смещения системы", "t, пс", "R², нм²");
        SetUpChart(ChartAcfSpeed, "График автокорреляционной функции скорости", "t, пс", "Z(t)");

        // Инициализация сцены для визуализации.
        _scene = new SceneManager { Viewport3D = Viewport };
    }

    #region ---СОБЫТИЯ СОЗДАНИЯ МОДЕЛИ---

    /// <summary>
    /// Событие создание модели.
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnCreateModel(object sender, RoutedEventArgs e)
    {
        if (_timer.Enabled)
            _timer.Stop();

        var size = NudSize.Value;
        var k = NudDisplacement.Value;

        BtnCreateModel.IsEnabled = false;
        ProgressBar.Value = 0;
        RtbOutputInfo.Document.Blocks.Clear();

        // Инициализация массивов энергий системы.
        _keValues = new List<double>();
        _peValues = new List<double>();
        _feValues = new List<double>();

        // Очистка графиков.
        ChartEnergy.Plot.Clear();
        ChartEnergy.Refresh();
        ChartMsd.Plot.Clear();
        ChartMsd.Refresh();
        ChartRadDist.Plot.Clear();
        ChartRadDist.Refresh();
        ChartAcfSpeed.Plot.Clear();
        ChartAcfSpeed.Refresh();

        _bgWorkerCreateModel.RunWorkerAsync(new FindPrimesInput(new object[] { size, k}));
    }

    /// <summary>
    /// DO WORK CREATE MODEL
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    /// <exception cref="ArgumentException"></exception>
    private void OnBackgroundWorkerDoWorkCreateModel(object sender, DoWorkEventArgs e)
    {
        var input = (FindPrimesInput)(e.Argument ?? throw new ArgumentException("Отсутствуют аргументы"));
        var size = (int)input.Args[0];
        var fraction = 0.5;
        var k = (double)input.Args[1];
        _bgWorkerCreateModel.ReportProgress(250);

        // Инициализация системы.
        _atomic = new AtomicModel(size, fraction);
        _initStep = _atomic.CurrentStep;
        _bgWorkerCreateModel.ReportProgress(500);

        // Применение случайного смещения для атомов.
        if (_isDisplacement)
        {
            _atomic.AtomsDisplacement(k);
            _bgWorkerCreateModel.ReportProgress(750);
        }

        // Вычисление начальных параметров системы.
        _atomic.InitCalculation();
        _bgWorkerCreateModel.ReportProgress(1000);

        // Начальное запоминание энергии системы.
        _keValues.Add(_atomic.Ke);
        _peValues.Add(_atomic.Pe);
        _feValues.Add(_atomic.Fe);
    }

    /// <summary>
    /// RUN WORKER COMPLETED CREATE MODEL
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnBackgroundWorkerRunWorkerCompletedCreateModel(object sender, RunWorkerCompletedEventArgs e)
    {
        if (e.Error != null)
            MessageBox.Show(e.Error.Message, "Произошла ошибка");
        else
        {
            _isNewSystem = true;

            // Запоминание позиции атомов на 0-ом шаге.
            _positionsAtomsList = new List<List<XYZ>> { _atomic.GetPosAtoms() };

            GC.Collect();
            // Отрисовка атомов на сцене.
            _scene.CreateScene(_positionsAtomsList.First(), _atomic.BoxSize, _atomic.GetSigma() / 2d);
            // Обнуление и блокировка слайдера.
            SliderTimeStep.Value = 0;
            SliderTimeStep.IsEnabled = false;

            // Вывод начальной информации.
            RtbOutputInfo.AppendText(InitInfoSystem());

            // Настройка и отрисовка графика радиального распределения.
            var rd = _atomic.GetRadialDistribution();
            _yMaxRb = rd.Max(p => p.Y);
            ChartRadDist.Plot.Clear();
            ChartRadDist.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(), Color.Blue, "Радиальное распределение");
            ChartRadDist.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.LatticeGeSn * 1e9 * 0.726, yMin: 0, yMax: _yMaxRb);
            ChartRadDist.Plot.Legend(location: Alignment.UpperRight);
            ChartRadDist.Refresh();

            BtnCreateModel.IsEnabled = true;
            BtnStartCalculation.IsEnabled = true;
            BtnCancelCalculation.IsEnabled = false;
            BtnToBegin.IsEnabled = false;
            BtnStepBack.IsEnabled = false;
            BtnPlayTimer.IsEnabled = false;
            BtnPauseTimer.IsEnabled = false;
            BtnStepForward.IsEnabled = false;
            BtnToEnd.IsEnabled = false;
            BtnFaster.IsEnabled = false;
            BtnSlower.IsEnabled = false;
        }
    }

    /// <summary>
    /// PROGRESS CHANGED CREATE MODEL
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnBackgroundWorkerProgressChangedCreateModel(object sender, ProgressChangedEventArgs e)
    {
        ProgressBar.Value = e.ProgressPercentage;
    }

    #endregion

    #region ---СОБЫТИЯ ЗАПУСКА МОДЕЛИРОВАНИЯ---

    /// <summary>
    /// Событие запуска/возобновления вычислений.
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnStartCalculation(object sender, RoutedEventArgs e)
    {
        var countStep = NudCountStep.Value;
        var snapshotStep = NudSnapshotStep.Value;
        var stepRt = NudStepMsd.Value;
        var T = NudTemperature.Value;
        var stepNorm = NudStepNorm.Value;

        _atomic.dt = (NudTimeStep.Value ?? 0.01) * 1e-12;
        _atomic.CountNumberAcf = (NudCountNumberAcf.Value ?? 150) + 1;
        _atomic.CountRepeatAcf = NudCountRepeatAcf.Value ?? 5;
        _atomic.StepRepeatAcf = NudStepRepeatAcf.Value ?? 10;

        BtnStartCalculation.IsEnabled = false;
        BtnCancelCalculation.IsEnabled = true;

        // Инициализация массива среднего квадрата смещения.
        _msdPoints = new List<PointD> { new(0, 0) };
        _averT = 0;

        // Очистка графиков.
        ChartEnergy.Plot.Clear();
        ChartEnergy.Plot.SetAxisLimits(_initStep * _atomic.dt * 1e12, (_initStep + countStep - 1) * _atomic.dt * 1e12);
        ChartEnergy.Refresh();
        ChartMsd.Plot.Clear();
        ChartMsd.Refresh();
        ChartAcfSpeed.Plot.Clear();
        ChartAcfSpeed.Refresh();

        // Сброс ProgressBar.
        ProgressBar.Value = 0;

        // Вывод начальной информации.
        RtbOutputInfo.AppendText(_isNormSpeeds ? "\n\nЗапуск перенормировки скоростей...\n" : "\n\nЗапуск моделирования...\n");
        RtbOutputInfo.AppendText($"Количество временных шагов: {countStep}\n" + (_isNormSpeeds ? $"Шаг перенормировки: {stepNorm}\n\n" : "\n"));
        RtbOutputInfo.AppendText(TableHeader());
        RtbOutputInfo.AppendText(TableData(_initStep - 1, 1));
        RtbOutputInfo.ScrollToEnd();

        // Запуск расчётов.
        _bgWorkerCalculation.RunWorkerAsync(new FindPrimesInput(new object[] { countStep, snapshotStep, stepRt, T, stepNorm }));
    }

    /// <summary>
    /// DO WORK CALCULATION
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    /// <exception cref="Exception"></exception>
    /// <exception cref="ArgumentException"></exception>
    private void OnBackgroundWorkerDoWorkCalculation(object sender, DoWorkEventArgs e)
    {
        if (_atomic == null)
            throw new NullReferenceException();

        var input = (FindPrimesInput)(e.Argument ?? throw new ArgumentException("Отсутствуют аргументы"));
        var countStep = (int)input.Args[0];
        var snapshotStep = (int)input.Args[1];
        var stepRt = (int)input.Args[2];
        var T = (int)input.Args[3];
        var stepNorm = (int)input.Args[4];

        // Начальная перенормировка скоростей, если она включено.
        if (_isNormSpeeds)
        {
            _atomic.InitVelocityNormalization(T);
            _atomic.PulseZeroing();
        }

        _atomic.CurrentStep = 1;

        // Запуск моделирования.
        for (var i = _initStep; i - _initStep < countStep; i++)
        {
            // Отслеживание отмены моделирования.
            if (_bgWorkerCalculation.CancellationPending)
            {
                e.Cancel = true;
                _initStep += _atomic.CurrentStep - 1;
                return;
            }

            // Расчёт шага методом Верле.
            _atomic.Verlet();

            // Проведение перенормировки скоростей, если она включено.
            if (_isNormSpeeds && i % stepNorm == 0)
                _atomic.VelocityNormalization(T);

            _keValues.Add(_atomic.Ke);
            _peValues.Add(_atomic.Pe);
            _feValues.Add(_atomic.Fe);
            _averT += _atomic.T;
            _averP += _atomic.P1;
            _positionsAtomsList.Add(_atomic.GetPosAtoms());

            // Вывод информации в UI.
            if ((_isSnapshot && i % snapshotStep == 0) || i == _initStep + countStep - 1)
                Application.Current.Dispatcher.Invoke(DispatcherPriority.Send, () =>
                {
                    RtbOutputInfo.AppendText(TableData(i, _isSnapshot ? snapshotStep : countStep));
                    RtbOutputInfo.ScrollToEnd();
                    _atomic.Flux = XYZ.Zero;

                    // Настройка и отрисовка графика радиального распределения.
                    var rd = _atomic.GetRadialDistribution();
                    ChartRadDist.Plot.Clear();
                    ChartRadDist.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(),
                        color: Color.Blue, label: "Радиальное распределение");
                    ChartRadDist.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.LatticeGeSn * 1e9 * 0.726, yMin: 0, yMax: _yMaxRb);
                    ChartRadDist.Plot.Legend(location: Alignment.UpperRight);
                    ChartRadDist.Refresh();
                });

            // Расчёт среднего квадрата смещения.
            if (i % stepRt == 0 || i == _initStep + countStep - 1)
                _msdPoints.Add(new PointD((i - _initStep + 1) * _atomic.dt, _atomic.GetMsd()));

            // Обновление ProgressBar.
            _bgWorkerCalculation.ReportProgress((int)((double)(i - _initStep) / countStep * 1000d) + 1);
        }

        _averT /= countStep;
        _averP /= countStep;
        _initStep += _atomic.CurrentStep - 1;
    }

    /// <summary>
    /// RUN WORKER COMPLETED CALCULATION
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnBackgroundWorkerRunWorkerCompletedCalculation(object sender, RunWorkerCompletedEventArgs e)
    {
        if (e.Cancelled)
            MessageBox.Show("Моделирование отменено");
        if (e.Error != null)
            MessageBox.Show(e.Error.Message, "Произошла ошибка");
        if (_atomic == null)
            throw new NullReferenceException();

        // Отрисовка графика энергий системы.
        ChartEnergy.Plot.AddSignal(_keValues.ToArray(), 1 / (_atomic.dt * 1e12), Color.Red, "Кинетическая энергия");
        ChartEnergy.Plot.AddSignal(_peValues.ToArray(), 1 / (_atomic.dt * 1e12), Color.Blue, "Потенциальная энергия");
        ChartEnergy.Plot.AddSignal(_feValues.ToArray(), 1 / (_atomic.dt * 1e12), Color.Green, "Полная энергия");
        ChartEnergy.Plot.AddHorizontalLine(0, Color.FromArgb(120, Color.Black));
        ChartEnergy.Plot.AddVerticalLine(0, Color.FromArgb(200, Color.Black));
        ChartEnergy.Plot.Margins(x: 0.0, y: 0.6);
        ChartEnergy.Plot.Legend(location: Alignment.UpperRight);
        ChartEnergy.Refresh();

        // Отрисовка графика радиального распределения.
        var rd = _atomic.GetRadialDistribution();
        ChartRadDist.Plot.Clear();
        ChartRadDist.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(), Color.Blue, "Радиальное распределение");
        ChartRadDist.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.LatticeGeSn * 1e9 * 0.726, yMin: 0, yMax: rd.Max(p => p.Y) * 1.1);
        ChartRadDist.Plot.Legend(location: Alignment.UpperRight);
        ChartRadDist.Refresh();

        // Отрисовка графика среднего квадрата смещения распределения.
        if (_msdPoints.Count != 1)
        {
            ChartMsd.Plot.AddSignalXY(_msdPoints.Select(p => p.X * 1e12).ToArray(), _msdPoints.Select(p => p.Y * 1e18).ToArray(), Color.Indigo, "Средний квадрат смещения");
            ChartMsd.Plot.SetAxisLimits(xMin: 0, xMax: _msdPoints.Max(p => p.X * 1e12), yMin: 0, yMax: (_msdPoints.Max(p => p.Y * 1e18) < 1e-10 ? 0.1 : _msdPoints.Max(p => p.Y * 1e18)) * 1.5);
            ChartMsd.Plot.Legend(location: Alignment.UpperRight);
            ChartMsd.Refresh();
        }

        // Отрисовка графика АКФ скорости.
        var zt = _atomic.GetAcfs(out var norm);
        ChartAcfSpeed.Plot.AddSignal(zt, 1 / (_atomic.dt * 1e12), Color.Green, "Автокорреляционная функция скорости");
        ChartAcfSpeed.Plot.SetAxisLimits(xMin: 0, xMax: (zt.Length - 1) * _atomic.dt * 1e12, yMin: -1, yMax: 1);
        ChartAcfSpeed.Plot.AddHorizontalLine(0, Color.FromArgb(120, Color.Black));
        ChartAcfSpeed.Plot.AddVerticalLine(0, Color.FromArgb(200, Color.Black));
        ChartAcfSpeed.Plot.Legend(location: Alignment.UpperRight);
        ChartAcfSpeed.Refresh();

        // Вывод информации в Rtb.
        var d1 = double.Round(_atomic.GetSelfDiffCoefFromAcf(zt, norm) * 1e9, 5);
        var d2 = double.Round(_atomic.GetSelfDiffCoefFromMsd(_msdPoints, out _) * 1e9, 5);
        var d3 = double.Round(_atomic.GetSelfDiffCoefFromMsd(_msdPoints[1], _msdPoints[_msdPoints.Count - 1]) * 1e9, 5);
        RtbOutputInfo.AppendText($"\n{double.Round(_averT, 3)} К - средняя температура");
        RtbOutputInfo.AppendText($"\n{double.Round(_averP, 1)} Па - среднее давление");
        RtbOutputInfo.AppendText($"\nDₛ ≈ {d1}•10⁻⁵ см²/с - коэф. самодифузии (полученный через АКФ)");
        RtbOutputInfo.AppendText($"\nDₛ ≈ {d2}•10⁻⁵ см²/с - коэф. самодифузии (полученный через средний квадрат смещения (МНК)) ");
        RtbOutputInfo.AppendText($"\nDₛ ≈ {d3}•10⁻⁵ см²/с - коэф. самодифузии (полученный через средний квадрат смещения (грубо)) ");

        _isNewSystem = false;
        BtnStartCalculation.IsEnabled = true;
        BtnCancelCalculation.IsEnabled = false;
        SliderTimeStep.IsEnabled = true;
        SliderTimeStep.Maximum = _positionsAtomsList.Count - 1;
        BtnToBegin.IsEnabled = false;
        BtnStepBack.IsEnabled = false;
        BtnPlayTimer.IsEnabled = true;
        BtnPauseTimer.IsEnabled = false;
        BtnStepForward.IsEnabled = true;
        BtnToEnd.IsEnabled = true;
        BtnFaster.IsEnabled = false;
        BtnSlower.IsEnabled = false;
    }

    /// <summary>
    /// PROGRESS CHANGED CALCULATION
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnBackgroundWorkerProgressChangedCalculation(object sender, ProgressChangedEventArgs e)
    {
        ProgressBar.Value = e.ProgressPercentage;
    }

    /// <summary>
    /// Событие отмены вычислений.
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnCancelCalculation(object sender, RoutedEventArgs e)
    {
        if (_bgWorkerCalculation.IsBusy)
            _bgWorkerCalculation.CancelAsync();
    }

    #endregion
}