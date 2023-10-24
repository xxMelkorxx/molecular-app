﻿using ScottPlot;
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
    private Dictionary<string, object> _params;
    private List<List<XYZ>> _positionsAtomsList;
    private List<PointD> _msdPoints;
    private double _averT, _averP;
    private bool _isDisplacement, _isSnapshot, _isNormSpeeds, _isNewSystem;
    private double _yMaxRb;
    private int _initStep;

    public MainWindow()
    {
        InitializeComponent();

        _bgWorkerCreateModel = (BackgroundWorker)FindResource("BackgroundWorkerCreateModel");
        _bgWorkerCalculation = (BackgroundWorker)FindResource("BackgroundWorkerCalculation");

        // Инициализация таймера для визуализации.
        _timer = new() { Interval = 30 };
        _timer.Tick += OnTickTimer;
    }

    private void OnLoadedMainWindow(object sender, RoutedEventArgs e)
    {
        // Настройка графиков.
        SetUpChart(Chart1, "Графики энергий системы", "t, пс", "E, эВ");
        SetUpChart(Chart3, "График радиального распределения системы", "r, нм", "g(r)");
        SetUpChart(Chart2, "График  среднего квадрата смещения системы", "t, пс", "R², нм²");
        SetUpChart(Chart4, "График автокорреляционной функции скорости", "t, пс", "Z(t)");

        _params = new Dictionary<string, object>(); // Инициализация словаря с параметрами.
        _scene = new SceneManager { Viewport3D = Viewport }; // Инициализация сцены для визуализации.
    }

    #region ---СОБЫТИЯ СОЗДАНИЯ МОДЕЛИ---

    /// <summary>
    /// Событие создание модели.
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    private void OnCreateModel(object sender, RoutedEventArgs e)
    {
        if (_timer.Enabled) // Остановка таймера.
            _timer.Stop();

        _params["size"] = NudSize.Value; // размер расчётной ячейки.
        _params["displacement"] = NudDisplacement.Value; // коэффициент начального смещения.
        _params["fraction"] = 0.5;
        // Инициализация массивов энергий системы.
        _params["ke"] = new List<double>();
        _params["pe"] = new List<double>();
        _params["fe"] = new List<double>();

        BtnCreateModel.IsEnabled = false;
        ProgressBar.Value = 0;
        RtbOutputInfo.Document.Blocks.Clear();

        // Очистка графиков.
        Chart1.Plot.Clear();
        Chart1.Refresh();
        Chart2.Plot.Clear();
        Chart2.Refresh();
        Chart3.Plot.Clear();
        Chart3.Refresh();
        Chart4.Plot.Clear();
        Chart4.Refresh();

        _bgWorkerCreateModel.RunWorkerAsync(); // Запуск создания модели.
    }

    /// <summary>
    /// DO WORK CREATE MODEL
    /// </summary>
    /// <param name="sender"></param>
    /// <param name="e"></param>
    /// <exception cref="ArgumentException"></exception>
    private void OnBackgroundWorkerDoWorkCreateModel(object sender, DoWorkEventArgs e)
    {
        // Инициализация системы.
        _atomic = new AtomicModel((int)_params["size"], (double)_params["fraction"]);
        _initStep = _atomic.CurrentStep;
        // Применение случайного смещения для атомов.
        if (_isDisplacement)
            _atomic.AtomsDisplacement((double)_params["displacement"]);
        // Вычисление начальных характеристик системы.
        _atomic.InitCalculation();
        ((List<double>)_params["ke"]).Add(_atomic.Ke);
        ((List<double>)_params["pe"]).Add(_atomic.Pe);
        ((List<double>)_params["fe"]).Add(_atomic.Fe);
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

            // Отрисовка атомов на сцене.
            _scene.CreateScene(_positionsAtomsList.First(), _atomic.BoxSize, _atomic.GetSigma() / 10d);

            // Вывод начальной информации.
            RtbOutputInfo.AppendText(InitInfoSystem());

            // Настройка и отрисовка графика радиального распределения.
            var rd = _atomic.GetRadialDistribution();
            _yMaxRb = rd.Max(p => p.Y);
            Chart3.Plot.Clear();
            Chart3.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(), Color.Blue, "Радиальное распределение");
            Chart3.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.LatticeGeSn * 1e9 * 0.726, yMin: 0, yMax: _yMaxRb);
            Chart3.Plot.Legend(location: Alignment.UpperRight);
            Chart3.Refresh();

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

            SliderTimeStep.Value = 0;
            SliderTimeStep.IsEnabled = false;

            GC.Collect();
        }
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
        BtnStartCalculation.IsEnabled = false;
        BtnCancelCalculation.IsEnabled = true;

        _params["countStep"] = NudCountStep.Value;
        _params["snapshotStep"] = NudSnapshotStep.Value;
        _params["stepRt"] = NudStepMsd.Value;
        _params["T"] = NudTemperature.Value;
        _params["stepNorm"] = NudStepNorm.Value;

        _atomic.dt = (NudTimeStep.Value ?? 0.01) * 1e-12;
        _atomic.CountNumberAcf = (NudCountNumberAcf.Value ?? 150) + 1;
        _atomic.CountRepeatAcf = NudCountRepeatAcf.Value ?? 5;
        _atomic.StepRepeatAcf = NudStepRepeatAcf.Value ?? 10;

        // Инициализация массива среднего квадрата смещения.
        _msdPoints = new List<PointD> { new(0, 0) };
        _averT = 0;

        // Очистка графиков.
        Chart1.Plot.Clear();
        Chart1.Plot.SetAxisLimits(
            xMin: _initStep * _atomic.dt * 1e12,
            xMax: (_initStep + (int)_params["countStep"] - 1) * _atomic.dt * 1e12);
        Chart1.Refresh();
        Chart2.Plot.Clear();
        Chart2.Refresh();
        Chart4.Plot.Clear();
        Chart4.Refresh();

        // Сброс ProgressBar.
        ProgressBar.Value = 0;
        ProgressBar.Maximum = (int)_params["countStep"];

        // Вывод начальной информации.
        RtbOutputInfo.AppendText(_isNormSpeeds ? "\n\nЗапуск перенормировки скоростей...\n" : "\n\nЗапуск моделирования...\n");
        RtbOutputInfo.AppendText($"Количество временных шагов: {(int)_params["countStep"]}\n" + (_isNormSpeeds ? $"Шаг перенормировки: {(int)_params["stepNorm"]}\n\n" : "\n"));
        RtbOutputInfo.AppendText(TableHeader());
        RtbOutputInfo.AppendText(TableData(_initStep - 1, 1));
        RtbOutputInfo.ScrollToEnd();

        // Запуск расчётов.
        _bgWorkerCalculation.RunWorkerAsync();
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

        var countStep = (int)_params["countStep"];
        var snapshotStep = (int)_params["snapshotStep"];

        _atomic.CurrentStep = 1;
        if (_isNormSpeeds) // Начальная перенормировка скоростей, если она включено.
        {
            _atomic.InitVelocityNormalization((double)_params["T"]);
            _atomic.PulseZeroing();
        }

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

            _atomic.Verlet(); // Расчёт шага методом Верле.

            if (_isNormSpeeds && i % (int)_params["stepNorm"] == 0) // Проведение перенормировки скоростей, если она включено.
                _atomic.VelocityNormalization((double)_params["T"]);

            ((List<double>)_params["ke"]).Add(_atomic.Ke);
            ((List<double>)_params["pe"]).Add(_atomic.Pe);
            ((List<double>)_params["fe"]).Add(_atomic.Fe);
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
                    Chart3.Plot.Clear();
                    Chart3.Plot.AddSignalXY(
                        rd.Select(p => p.X * 1e9).ToArray(),
                        rd.Select(p => p.Y).ToArray(),
                        color: Color.Blue, label: "Радиальное распределение"
                    );
                    Chart3.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.LatticeGeSn * 1e9 * 0.726, yMin: 0, yMax: _yMaxRb);
                    Chart3.Plot.Legend(location: Alignment.UpperRight);
                    Chart3.Refresh();
                });

            // Расчёт среднего квадрата смещения.
            if (i % (int)_params["stepRt"] == 0 || i == _initStep + countStep - 1)
                _msdPoints.Add(new PointD((i - _initStep + 1) * _atomic.dt, _atomic.GetMsd()));

            // Обновление ProgressBar.
            _bgWorkerCalculation.ReportProgress(i - _initStep);
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
        Chart1.Plot.AddSignal(((List<double>)_params["ke"]).ToArray(), 1 / (_atomic.dt * 1e12), Color.Red, "Кинетическая энергия");
        Chart1.Plot.AddSignal(((List<double>)_params["pe"]).ToArray(), 1 / (_atomic.dt * 1e12), Color.Blue, "Потенциальная энергия");
        Chart1.Plot.AddSignal(((List<double>)_params["fe"]).ToArray(), 1 / (_atomic.dt * 1e12), Color.Green, "Полная энергия");
        Chart1.Plot.AddHorizontalLine(0, Color.FromArgb(120, Color.Black));
        Chart1.Plot.AddVerticalLine(0, Color.FromArgb(200, Color.Black));
        Chart1.Plot.Margins(x: 0.0, y: 0.6);
        Chart1.Plot.Legend(location: Alignment.UpperRight);
        Chart1.Refresh();

        // Отрисовка графика радиального распределения.
        // var rd = _atomic.GetRadialDistribution();
        // Chart3.Plot.Clear();
        // Chart3.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(), Color.Blue, "Радиальное распределение");
        // Chart3.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.LatticeGeSn * 1e9 * 0.726, yMin: 0, yMax: rd.Max(p => p.Y) * 1.1);
        // Chart3.Plot.Legend(location: Alignment.UpperRight);
        // Chart3.Refresh();

        // Отрисовка графика среднего квадрата смещения распределения.
        // if (_msdPoints.Count != 1)
        // {
        //     Chart2.Plot.AddSignalXY(_msdPoints.Select(p => p.X * 1e12).ToArray(), _msdPoints.Select(p => p.Y * 1e18).ToArray(), Color.Indigo, "Средний квадрат смещения");
        //     Chart2.Plot.SetAxisLimits(xMin: 0, xMax: _msdPoints.Max(p => p.X * 1e12), yMin: 0, yMax: (_msdPoints.Max(p => p.Y * 1e18) < 1e-10 ? 0.1 : _msdPoints.Max(p => p.Y * 1e18)) * 1.5);
        //     Chart2.Plot.Legend(location: Alignment.UpperRight);
        //     Chart2.Refresh();
        // }

        // Отрисовка графика АКФ скорости.
        // var zt = _atomic.GetAcfs(out var norm);
        // Chart4.Plot.AddSignal(zt, 1 / (_atomic.dt * 1e12), Color.Green, "Автокорреляционная функция скорости");
        // Chart4.Plot.SetAxisLimits(xMin: 0, xMax: (zt.Length - 1) * _atomic.dt * 1e12, yMin: -1, yMax: 1);
        // Chart4.Plot.AddHorizontalLine(0, Color.FromArgb(120, Color.Black));
        // Chart4.Plot.AddVerticalLine(0, Color.FromArgb(200, Color.Black));
        // Chart4.Plot.Legend(location: Alignment.UpperRight);
        // Chart4.Refresh();

        // Вывод информации в Rtb.
        // var d1 = double.Round(_atomic.GetSelfDiffCoefFromAcf(zt, norm) * 1e9, 5);
        // var d2 = double.Round(_atomic.GetSelfDiffCoefFromMsd(_msdPoints, out _) * 1e9, 5);
        // var d3 = double.Round(_atomic.GetSelfDiffCoefFromMsd(_msdPoints[1], _msdPoints[_msdPoints.Count - 1]) * 1e9, 5);
        // RtbOutputInfo.AppendText($"\n{double.Round(_averT, 3)} К - средняя температура");
        // RtbOutputInfo.AppendText($"\n{double.Round(_averP, 1)} Па - среднее давление");
        // RtbOutputInfo.AppendText($"\nDₛ ≈ {d1}•10⁻⁵ см²/с - коэф. самодифузии (полученный через АКФ)");
        // RtbOutputInfo.AppendText($"\nDₛ ≈ {d2}•10⁻⁵ см²/с - коэф. самодифузии (полученный через средний квадрат смещения (МНК)) ");
        // RtbOutputInfo.AppendText($"\nDₛ ≈ {d3}•10⁻⁵ см²/с - коэф. самодифузии (полученный через средний квадрат смещения (грубо)) ");

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