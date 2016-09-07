package org.processmining.antialignments.ilp.antialignment;

import info.clearthought.layout.TableLayout;
import info.clearthought.layout.TableLayoutConstants;

import java.awt.Dimension;

import javax.swing.JComponent;

import com.fluxicon.slickerbox.components.NiceDoubleSlider;
import com.fluxicon.slickerbox.components.NiceIntegerSlider;
import com.fluxicon.slickerbox.components.NiceSlider.Orientation;
import com.fluxicon.slickerbox.factory.SlickerFactory;

public class AntiAlignmentParameterUI extends JComponent {

	private static final long serialVersionUID = -6264123593622120018L;
	private final NiceIntegerSlider cutOffLengthSlider;
	private final NiceDoubleSlider maxFactorSlider;

	public AntiAlignmentParameterUI() {
		TableLayout tl = new TableLayout(new double[][] { { TableLayoutConstants.FILL }, { 80, 40, 40 } });

		setLayout(tl);

		SlickerFactory slickerFactoryInstance = SlickerFactory.instance();
		setTitle(
				slickerFactoryInstance,
				"<html><h1>Set parameters</h1><p>Use the sliders to set the prefix length and the maximum anti-alignment length factor.</p></html>");

		// max modelMoves
		cutOffLengthSlider = slickerFactoryInstance.createNiceIntegerSlider("<html><h4>Prefix length</h4></html>", 0,
				100, 5, Orientation.HORIZONTAL);
		cutOffLengthSlider.setPreferredSize(new Dimension(700, 20));
		cutOffLengthSlider.setMaximumSize(new Dimension(700, 20));
		cutOffLengthSlider.setMinimumSize(new Dimension(700, 20));
		add(cutOffLengthSlider, "0, 1, c, t");

		maxFactorSlider = slickerFactoryInstance.createNiceDoubleSlider("<html><h4>Maximum length factor</h4></html>",
				1.0, 5.0, 1.0, Orientation.HORIZONTAL);
		maxFactorSlider.setPreferredSize(new Dimension(700, 20));
		maxFactorSlider.setMaximumSize(new Dimension(700, 20));
		maxFactorSlider.setMinimumSize(new Dimension(700, 20));
		add(maxFactorSlider, "0, 2, c, t");

	}

	protected void setTitle(SlickerFactory slickerFactoryInstance, String title) {
		add(slickerFactoryInstance.createLabel(title), "0, 0, l, t");
	}

	public double getMaxFactor() {
		return maxFactorSlider.getValue();
	}

	public int getCutoffLength() {
		return cutOffLengthSlider.getValue();
	}

	public AntiAlignmentParameters getParameters() {
		return new AntiAlignmentParameters(getCutoffLength(), getMaxFactor());
	}
}
