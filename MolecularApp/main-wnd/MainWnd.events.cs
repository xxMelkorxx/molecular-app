using System.Windows;

namespace MolecularApp;

public partial class MainWindow
{
    private void OnCheckedDisplacement(object sender, RoutedEventArgs e)
    {
        _isDisplacement = true;
        NudDisplacement.IsEnabled = true;
    }
    
    private void OnUncheckedDisplacement(object sender, RoutedEventArgs e)
    {
        _isDisplacement = false;
        NudDisplacement.IsEnabled = false;
    }

    private void OnCheckedIsSnapshot(object sender, RoutedEventArgs e)
    {
        _isSnapshot = true;
        NudSnapshotStep.IsEnabled = true;
    }

    private void OnUncheckedIsSnapshot(object sender, RoutedEventArgs e)
    {
        _isSnapshot = false;
        NudSnapshotStep.IsEnabled = false;
    }

    private void OnCheckedIsRenormSpeeds(object sender, RoutedEventArgs e)
    {
        _isNormSpeeds = true;
        NudTemperature.IsEnabled = true;
        NudStepNorm.IsEnabled = true;
    }

    private void OnUncheckedIsRenormSpeeds(object sender, RoutedEventArgs e)
    {
        _isNormSpeeds = false;
        NudTemperature.IsEnabled = false;
        NudStepNorm.IsEnabled = false;
    }

    private void OnValueChangedNudCountNumberAcf(object sender, RoutedPropertyChangedEventArgs<object> e)
    {
        if (NudCountNumberAcf.Value > NudCountStep.Value) NudCountNumberAcf.Value = NudCountStep.Value;
    }
    
    private void OnValueChangedNudSnapshotStep(object sender, RoutedPropertyChangedEventArgs<object> e)
    {
        if (NudSnapshotStep.Value > NudCountStep.Value) NudSnapshotStep.Value = NudCountStep.Value;
    }

    private void OnValueChangedNudCountStep(object sender, RoutedPropertyChangedEventArgs<object> e)
    {
        if (NudCountNumberAcf is not null && NudCountNumberAcf.Value > NudCountStep.Value) NudCountNumberAcf.Value = NudCountStep.Value;
        if (NudSnapshotStep is not null && NudSnapshotStep.Value > NudCountStep.Value) NudSnapshotStep.Value = NudCountStep.Value;
    }
}