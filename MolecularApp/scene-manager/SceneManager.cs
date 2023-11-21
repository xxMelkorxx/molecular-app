using HelixToolkit.Wpf;
using System.Collections.Generic;
using System.Windows.Media;
using System.Windows.Media.Media3D;

namespace MolecularApp.scene_manager;

public class SceneManager
{
    /// <summary>
    /// Элемент, в которой создаётся сцена.
    /// </summary>
    public HelixViewport3D Viewport3D { get; set; }

    private List<XYZ> _initPosAtoms;

    /// <summary>
    /// Отрисовка сцены.
    /// </summary>
    /// <param name="positons">Список координат атомов</param>
    /// <param name="l">Размер расчётной ячейки (м)</param>
    /// <param name="radius">Радиус атомов (м)</param>
    public void CreateScene(List<XYZ> positons, double l, double radius)
    {
        var posCamera = new Point3D(-l * 3 * 1e9, l * 1e9, -l * 2 * 1e9);
        var rotateCenter = new Point3D(0, 0, 0);
        var dirCamera = rotateCenter - posCamera;

        // Добавление камеры на сцену.
        PerspectiveCamera camera = new()
        {
            Position = posCamera,
            LookDirection = dirCamera,
        };
        Viewport3D.Camera = camera;
        Viewport3D.Items.Clear();

        // Добавление источника света на сцену.
        Viewport3D.Children.Add(new DefaultLights());

        // Запоминание начального положения атомов.
        _initPosAtoms = new List<XYZ>();

        // Размещение атомов на сцену.
        positons.ForEach(pos =>
        {
            _initPosAtoms.Add(pos - l / 2);
            SphereVisual3D sphere = new()
            {
                Center = new Point3D((pos.X - l / 2) * 1e9, (pos.Y - l / 2) * 1e9, (pos.Z - l / 2) * 1e9),
                Radius = radius * 1e9,
                Material = Materials.Blue
            };
            Viewport3D.Items.Add(sphere);
        });

        // Добавление полупрозрачной коробки, визуалирующая размеры расчётной ячейки.
        BoxVisual3D box = new()
        {
            Length = l * 1e9,
            Width = l * 1e9,
            Height = l * 1e9,
            Center = rotateCenter,
            Material = new DiffuseMaterial(new SolidColorBrush(Color.FromArgb(100, 100, 100, 100)))
        };
        Viewport3D.Items.Add(box);
    }

    /// <summary>
    /// Обновление положения атомов.
    /// </summary>
    /// <param name="positons">Список координат атомов.</param>
    /// <param name="l">Размеры расчётной ячейки.</param>
    public void UpdatePositionsAtoms(List<XYZ> positons, double l)
    {
        for (var i = 0; i < positons.Count; i++)
            ((SphereVisual3D)Viewport3D.Items[i]).Transform = new TranslateTransform3D(
                (positons[i].X - l / 2 - _initPosAtoms[i].X) * 1e9,
                (positons[i].Y - l / 2 - _initPosAtoms[i].Y) * 1e9,
                (positons[i].Z - l / 2 - _initPosAtoms[i].Z) * 1e9
            );
        ((BoxVisual3D)Viewport3D.Items[positons.Count]).Length = l * 1e9;
        ((BoxVisual3D)Viewport3D.Items[positons.Count]).Width = l * 1e9;
        ((BoxVisual3D)Viewport3D.Items[positons.Count]).Height = l * 1e9;
    }
}