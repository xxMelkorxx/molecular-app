using System;
using System.Windows;

namespace MolecularApp;

public partial class MainWindow
{
    private void OnValueChangedSliderTimeStep(object sender, RoutedPropertyChangedEventArgs<double> e)
    {
        if (_atomic == null)
            throw new NullReferenceException();

        if ((int)SliderTimeStep.Value == (int)SliderTimeStep.Minimum)
        {
            BtnToBegin.IsEnabled = false;
            BtnStepBack.IsEnabled = false;
            BtnToEnd.IsEnabled = true;
            BtnStepForward.IsEnabled = true;
        }
        else if ((int)SliderTimeStep.Value == (int)SliderTimeStep.Maximum)
        {
            BtnToBegin.IsEnabled = true;
            BtnStepBack.IsEnabled = true;
            BtnToEnd.IsEnabled = false;
            BtnStepForward.IsEnabled = false;
        }
        else
        {
            BtnToBegin.IsEnabled = true;
            BtnStepBack.IsEnabled = true;
            BtnToEnd.IsEnabled = true;
            BtnStepForward.IsEnabled = true;
        }

        if (_isNewSystem)
            _scene.CreateScene(_positionsAtomsList[0], _boxSizeValues[0], _atomic.GetRadiusAtom() / 2);
        else
            _scene.UpdatePositionsAtoms(_positionsAtomsList[(int)SliderTimeStep.Value], _boxSizeValues[(int)SliderTimeStep.Value]);
    }

    private void OnClickBtnToBegin(object sender, RoutedEventArgs e)
    {
        SliderTimeStep.Value = SliderTimeStep.Minimum;
    }

    private void OnClickBtnStepBack(object sender, RoutedEventArgs e)
    {
        SliderTimeStep.Value--;
    }

    private void OnClickBtnPlayTimer(object sender, RoutedEventArgs e)
    {
        _timer.Start();

        BtnPlayTimer.Visibility = Visibility.Collapsed;
        BtnPauseTimer.Visibility = Visibility.Visible;
        BtnFaster.IsEnabled = _timer.Interval != 1;
        BtnSlower.IsEnabled = true;
    }

    private void OnClickBtnPauseTimer(object sender, RoutedEventArgs e)
    {
        _timer.Stop();

        BtnPlayTimer.Visibility = Visibility.Visible;
        BtnPauseTimer.Visibility = Visibility.Collapsed;
        BtnFaster.IsEnabled = false;
        BtnSlower.IsEnabled = false;
    }

    private void OnClickBtnStepForward(object sender, RoutedEventArgs e)
    {
        SliderTimeStep.Value++;
    }

    private void OnClickBtnToEnd(object sender, RoutedEventArgs e)
    {
        SliderTimeStep.Value = SliderTimeStep.Maximum;
    }

    private void OnClickBtnSlower(object sender, RoutedEventArgs e)
    {
        _timer.Interval += 5;
        BtnFaster.IsEnabled = true;
    }

    private void OnClickBtnFaster(object sender, RoutedEventArgs e)
    {
        if (_timer.Interval - 5 <= 1)
        {
            _timer.Interval = 1;
            BtnFaster.IsEnabled = false;
        }
        else
            _timer.Interval -= 5;
    }

    private void OnTickTimer(object sender, EventArgs e)
    {
        if ((int)SliderTimeStep.Value == (int)SliderTimeStep.Maximum)
            SliderTimeStep.Value = SliderTimeStep.Minimum;
        else
            SliderTimeStep.Value++;
    }
}