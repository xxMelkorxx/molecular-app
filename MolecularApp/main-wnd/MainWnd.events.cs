using System.Windows;
using System.Windows.Controls;

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

    private void OnValueChangedFirstFraction(object sender, RoutedPropertyChangedEventArgs<object> e)
    {
        if (e.OldValue != null)
            NudSecondFraction.Value = 1d - (double)e.NewValue;
    }

    private void OnValueChangedSecondFraction(object sender, RoutedPropertyChangedEventArgs<object> e)
    {
        if (e.OldValue != null)
            NudFirstFraction.Value = 1d - (double)e.NewValue;
    }

    private void OnSelectionChangedComboBoxAtomsType(object sender, SelectionChangedEventArgs e) { _atomType = (AtomType)ComboBoxAtomsType.SelectedIndex; }

    private void OnSelectionChangedComboBoxFirstAtomsType(object sender, SelectionChangedEventArgs e) { _firstAtom = (AtomType)ComboBoxFirstAtomsType.SelectedIndex; }

    private void OnSelectionChangedComboBoxSecondAtomsType(object sender, SelectionChangedEventArgs e) { _secondAtom = (AtomType)ComboBoxSecondAtomsType.SelectedIndex; }

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

    private void OnCheckedRadioButtonCrystal(object sender, RoutedEventArgs e)
    {
        if (CrystalParamsPanel is not null)
        {
            _isCrystal = true;
            CrystalParamsPanel.Visibility = Visibility.Visible;
            AlloyParamsPanel.Visibility = Visibility.Collapsed;
        }
    }

    private void OnUncheckedRadioButtonCrystal(object sender, RoutedEventArgs e)
    {
        if (AlloyParamsPanel is not null)
        {
            _isCrystal = false;
            CrystalParamsPanel.Visibility = Visibility.Collapsed;
            AlloyParamsPanel.Visibility = Visibility.Visible;
        }
    }
}