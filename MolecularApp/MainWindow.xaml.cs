using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Drawing;
using System.IO;
using System.Windows;
using System.Windows.Threading;
using ScottPlot;
using MolecularApp.atomic_model;
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
    private int _initStep, _iter;
    private AtomType _firstAtom, _secondAtom;
    private double _l0;

    private const string saveImagePath = "S:\\SerBor\\Научка\\ImageResults";

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
        SetUpChart(Chart2, "График радиального распределения системы", "r, нм", "g(r)");
        SetUpChart(Chart3, "График  среднего квадрата смещения системы", "t, пс", "R², нм²");
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
        _params["firstFraction"] = NudFirstFraction.Value;
        _params["secondFraction"] = NudSecondFraction.Value;
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
        Chart3.Plot.Clear();
        Chart3.Refresh();
        Chart2.Plot.Clear();
        Chart2.Refresh();
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
        _atomic = new AtomicModel(
            size: (int)_params["size"],
            firstTypeAtom: _firstAtom,
            fisrtFraction: (double)_params["firstFraction"],
            secondTypeAtom: _secondAtom,
            secondFraction: (double)_params["secondFraction"]
        );
        _initStep = _atomic.CurrentStep;
        _l0 = _atomic.SystemLattice;

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
            _scene.CreateScene(_positionsAtomsList.First(), _atomic.BoxSize, _atomic.GetRadiusAtom() / 2d);

            // Вывод начальной информации.
            RtbOutputInfo.AppendText(InitInfoSystem());

            // Настройка и отрисовка графика радиального распределения.
            var rd = _atomic.GetRadialDistribution();
            _yMaxRb = rd.Max(p => p.Y) * 1.5;
            Chart2.Plot.Clear();
            Chart2.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(), Color.Blue, "Радиальное распределение");
            Chart2.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.SystemLattice * 1e9 * 0.726, yMin: 0, yMax: _yMaxRb);
            Chart2.Plot.Legend(location: Alignment.UpperRight);
            Chart2.Refresh();

            BtnCreateModel.IsEnabled = true;
            BtnStartCalculation.IsEnabled = true;
            BtnCancelCalculation.IsEnabled = false;
            BtnToBegin.IsEnabled = false;
            BtnStepBack.IsEnabled = false;
            BtnPlayTimer.IsEnabled = false;
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
        _params["tempStep"] = NudTempStep.Value;
        _params["snapshotStep"] = NudSnapshotStep.Value;
        _params["stepRt"] = NudStepMsd.Value;
        _params["T"] = NudTemperature.Value;
        _params["stepNorm"] = NudStepNorm.Value;

        // Расширение расчётной ячейки от температуры.
        var coefUp = NudUpperBoxSizeСoef.Value ?? 1d;
        _atomic.SystemLattice = _l0 * coefUp;
        _atomic.dt = (NudTimeStep.Value ?? 0.01) * 1e-12;
        _atomic.CountNumberAcf = (NudCountNumberAcf.Value ?? 150) + 1;
        _atomic.CountRepeatAcf = NudCountRepeatAcf.Value ?? 5;
        _atomic.StepRepeatAcf = NudStepRepeatAcf.Value ?? 10;
        _atomic.PulseZeroing();

        // Инициализация массива среднего квадрата смещения.
        _msdPoints = new List<PointD> { new(0, 0) };
        _averT = 0;
        _averP = 0;
        _iter = 0;
        
        // Очистка графиков.
        Chart1.Plot.Clear();
        Chart1.Plot.SetAxisLimits(
            xMin: _initStep * _atomic.dt * 1e12,
            xMax: (_initStep + (int)_params["countStep"] - 1) * _atomic.dt * 1e12);
        Chart1.Refresh();
        Chart3.Plot.Clear();
        Chart3.Refresh();
        Chart4.Plot.Clear();
        Chart4.Refresh();

        // Сброс ProgressBar.
        ProgressBar.Value = 0;
        ProgressBar.Maximum = (int)_params["countStep"] - 1;

        // Вывод начальной информации.
        RtbOutputInfo.AppendText(_isNormSpeeds ? "\n\nЗапуск перенормировки скоростей...\n" : "\n\nЗапуск моделирования...\n");
        RtbOutputInfo.AppendText($"Количество временных шагов: {(int)_params["countStep"]}\n" + (_isNormSpeeds ? $"Шаг перенормировки: {(int)_params["stepNorm"]}\n\n" : "\n"));
        RtbOutputInfo.AppendText(TableHeader());
        RtbOutputInfo.AppendText(TableData(_initStep - 1, 1));
        RtbOutputInfo.ScrollToEnd();

        // Создание каталога, для сохранения данных.
        if (_isNewSystem)
        {
            var saveDirectory = new DirectoryInfo(
                $"{saveImagePath}\\Results_{DateTime.Now.ToString("ddmmyyyy_hhmmss")}_{_firstAtom}{_secondAtom}_{(int)(_atomic.SecondFraction * 100)}_{(int)(_atomic.FisrtFraction * 100)}_{_atomic.CountAtoms}"
            );
            if (!saveDirectory.Exists)
                saveDirectory.Create();
            saveDirectory.CreateSubdirectory("Energy");
            saveDirectory.CreateSubdirectory("Rad");
            saveDirectory.CreateSubdirectory("Msd");
            saveDirectory.CreateSubdirectory("Acf");
            _params["saveDirectory"] = saveDirectory.FullName;
        }
        
        // Запоминание позиций атомов.
        _positionsAtomsList = new List<List<XYZ>> { _atomic.GetPosAtoms() };
        GC.Collect();
        SliderTimeStep.Value = 0;
        
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
        var tempStep = (int)_params["tempStep"];
        var snapshotStep = (int)_params["snapshotStep"];

        _atomic.CurrentStep = 1;

        // Начальная перенормировка скоростей, если она включено.
        if (_isNormSpeeds)
            _atomic.InitVelocityNormalization((double)_params["T"]);

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
            if (_isNormSpeeds && i % (int)_params["stepNorm"] == 0)
                _atomic.VelocityNormalization((double)_params["T"]);

            ((List<double>)_params["ke"]).Add(_atomic.Ke);
            ((List<double>)_params["pe"]).Add(_atomic.Pe);
            ((List<double>)_params["fe"]).Add(_atomic.Fe);
            _averT += _atomic.T;
            _averP += _atomic.P1;
            _iter++;
            _positionsAtomsList.Add(_atomic.GetPosAtoms());

            // Вывод информации в UI.
            if ((_isSnapshot && i % snapshotStep == 0) || i == _initStep + countStep - 1)
                Application.Current.Dispatcher.Invoke(DispatcherPriority.Send, () =>
                {
                    RtbOutputInfo.AppendText(TableData(i, _isSnapshot ? snapshotStep : countStep));
                    RtbOutputInfo.ScrollToEnd();

                    // Настройка и отрисовка графика радиального распределения.
                    var rd = _atomic.GetRadialDistribution();
                    Chart2.Plot.Clear();
                    Chart2.Plot.AddSignalXY(
                        rd.Select(p => p.X * 1e9).ToArray(),
                        rd.Select(p => p.Y).ToArray(),
                        color: Color.Blue, label: "Радиальное распределение"
                    );
                    Chart2.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.SystemLattice * 1e9 * 0.726, yMin: 0, yMax: _yMaxRb);
                    Chart2.Plot.Legend(location: Alignment.UpperRight);
                    Chart2.Refresh();
                });

            if (_iter == tempStep)
            {
                Application.Current.Dispatcher.Invoke(DispatcherPriority.Send, () =>
                {
                    RtbOutputInfo.AppendText($"\n{(_averT / tempStep).ToString("F1")} К - средняя температура");
                    RtbOutputInfo.AppendText($"\n{(_averP / tempStep / 1e6).ToString("F1")} МПа - среднее давление (через вириал)");
                    RtbOutputInfo.AppendText($"\n{(_atomic.P2 / tempStep / 1e6).ToString("F1")} МПа - среднее давление\n\n");
                    RtbOutputInfo.ScrollToEnd();
                });
                _averT = 0;
                _averP = 0;
                _atomic.Flux = XYZ.Zero;
                _iter = 0;
            }

            // Расчёт среднего квадрата смещения.
            if (i % (int)_params["stepRt"] == 0 || i == _initStep + countStep - 1)
                _msdPoints.Add(new PointD((i - _initStep + 1) * _atomic.dt, _atomic.GetMsd()));

            // Обновление ProgressBar.
            _bgWorkerCalculation.ReportProgress(i - _initStep);
        }

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
        Chart1.Plot.SaveFig(
            $"{_params["saveDirectory"]}\\Energy\\Steps_{_initStep - 1}_T_{((double)_params["T"]).ToString("F1")}.png",
            width: 1500, height: 1200
        );

        // Отрисовка графика радиального распределения.
        var rd = _atomic.GetRadialDistribution();
        Chart2.Plot.Clear();
        Chart2.Plot.AddSignalXY(rd.Select(p => p.X * 1e9).ToArray(), rd.Select(p => p.Y).ToArray(), Color.Blue, "Радиальное распределение");
        Chart2.Plot.SetAxisLimits(xMin: 0, xMax: 5 * _atomic.SystemLattice * 1e9 * 0.726, yMin: 0, yMax: rd.Max(p => p.Y) * 1.1);
        Chart2.Plot.Legend(location: Alignment.UpperRight);
        Chart2.Refresh();
        Chart2.Plot.SaveFig(
            $"{_params["saveDirectory"]}\\Rad\\Steps_{_initStep - 1}_T_{((double)_params["T"]).ToString("F0")}.png",
            width: 1500, height: 1200
        );

        // Отрисовка графика среднего квадрата смещения распределения.
        if (_msdPoints.Count != 1)
        {
            Chart3.Plot.AddSignalXY(_msdPoints.Select(p => p.X * 1e12).ToArray(), _msdPoints.Select(p => p.Y * 1e18).ToArray(), Color.Indigo, "Средний квадрат смещения");
            Chart3.Plot.SetAxisLimits(xMin: 0, xMax: _msdPoints.Max(p => p.X * 1e12), yMin: 0, yMax: (_msdPoints.Max(p => p.Y * 1e18) < 1e-10 ? 0.1 : _msdPoints.Max(p => p.Y * 1e18)) * 1.5);
            Chart3.Plot.Legend(location: Alignment.UpperRight);
            Chart3.Refresh();
            Chart3.Plot.SaveFig(
                $"{_params["saveDirectory"]}\\Msd\\Steps_{_initStep - 1}_T_{((double)_params["T"]).ToString("F0")}.png",
                width: 1500, height: 1200
            );
        }

        // Отрисовка графика АКФ скорости.
        var zt = _atomic.GetAcfs(out var norm);
        Chart4.Plot.AddSignal(zt, 1 / (_atomic.dt * 1e12), Color.Green, "Автокорреляционная функция скорости");
        Chart4.Plot.SetAxisLimits(xMin: 0, xMax: (zt.Length - 1) * _atomic.dt * 1e12, yMin: -1, yMax: 1);
        Chart4.Plot.AddHorizontalLine(0, Color.FromArgb(120, Color.Black));
        Chart4.Plot.AddVerticalLine(0, Color.FromArgb(200, Color.Black));
        Chart4.Plot.Legend(location: Alignment.UpperRight);
        Chart4.Refresh();
        Chart4.Plot.SaveFig(
            $"{_params["saveDirectory"]}\\Acf\\Steps_{_initStep - 1}_T_{((double)_params["T"]).ToString("F0")}.png",
            width: 1500, height: 1200
        );

        // Вывод информации в Rtb.
        var d1 = double.Round(_atomic.GetSelfDiffCoefFromAcf(zt, norm) * 1e9, 5);
        var d2 = double.Round(_atomic.GetSelfDiffCoefFromMsd(_msdPoints, out _) * 1e9, 5);
        var d3 = double.Round(_atomic.GetSelfDiffCoefFromMsd(_msdPoints[1], _msdPoints[_msdPoints.Count - 1]) * 1e9, 5);
        RtbOutputInfo.AppendText($"Dₛ ≈ {d1}•10⁻⁵ см²/с - коэф. самодифузии (полученный через АКФ)\n");
        RtbOutputInfo.AppendText($"Dₛ ≈ {d2}•10⁻⁵ см²/с - коэф. самодифузии (полученный через средний квадрат смещения (МНК))\n");
        RtbOutputInfo.AppendText($"Dₛ ≈ {d3}•10⁻⁵ см²/с - коэф. самодифузии (полученный через средний квадрат смещения (грубо))\n");

        // Звуковое оповещение.
        AlarmBeep(500, 500, 1);

        _isNewSystem = false;
        BtnStartCalculation.IsEnabled = true;
        BtnCancelCalculation.IsEnabled = false;
        SliderTimeStep.IsEnabled = true;
        SliderTimeStep.Maximum = _positionsAtomsList.Count - 1;
        BtnToBegin.IsEnabled = false;
        BtnStepBack.IsEnabled = false;
        BtnPlayTimer.IsEnabled = true;
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
    private void OnBackgroundWorkerProgressChangedCalculation(object sender, ProgressChangedEventArgs e) { ProgressBar.Value = e.ProgressPercentage; }

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