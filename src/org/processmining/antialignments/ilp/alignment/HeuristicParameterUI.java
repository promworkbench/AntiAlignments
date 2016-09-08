package org.processmining.antialignments.ilp.alignment;

import info.clearthought.layout.TableLayoutConstants;

import java.awt.Dimension;
import java.util.Collection;
import java.util.Map;

import org.deckfour.xes.classification.XEventClass;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.plugins.petrinet.replayer.algorithms.costbasedcomplete.CostBasedCompleteUI;

import com.fluxicon.slickerbox.components.NiceDoubleSlider;
import com.fluxicon.slickerbox.components.NiceIntegerSlider;
import com.fluxicon.slickerbox.components.NiceSlider.Orientation;
import com.fluxicon.slickerbox.factory.SlickerFactory;

public class HeuristicParameterUI extends CostBasedCompleteUI {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1750372583381967369L;
	private static final int MODELMOVES = 5;
	private NiceIntegerSlider modelMoves;
	private NiceIntegerSlider backtrackMax;
	private NiceDoubleSlider backtrackThreshold;

	public HeuristicParameterUI(Collection<Transition> transCol, Collection<XEventClass> evClassCol, int modelSize) {
		this(transCol, evClassCol, null, null, null, modelSize);
	}

	public HeuristicParameterUI(Collection<Transition> transCol, Collection<XEventClass> evClassCol,
			Map<Transition, Integer> defMoveModelCost, Map<Transition, Integer> defSyncCost,
			Map<XEventClass, Integer> defMoveLogCost, int modelSize) {
		super(new double[][] {
				{ TableLayoutConstants.FILL },
				{ 80, 40, 40, 40, TableLayoutConstants.FILL, 35, TableLayoutConstants.FILL, 35,
						TableLayoutConstants.FILL, 35 } });

		SlickerFactory slickerFactoryInstance = SlickerFactory.instance();
		setTitle(
				slickerFactoryInstance,
				"<html><h1>Set parameters</h1><p>Double click costs on table to change their values. Use only non-negative integers.</p></html>");

		// max modelMoves
		modelMoves = slickerFactoryInstance.createNiceIntegerSlider(
				"<html><h4># Expected maximal model move sequence.</h4></html>", 0, 2 * modelSize,
				Math.min(modelSize, MODELMOVES), Orientation.HORIZONTAL);
		modelMoves.setPreferredSize(new Dimension(700, 20));
		modelMoves.setMaximumSize(new Dimension(700, 20));
		modelMoves.setMinimumSize(new Dimension(700, 20));
		add(modelMoves, "0, 1, c, t");

		// max backtrack
		backtrackMax = slickerFactoryInstance.createNiceIntegerSlider(
				"<html><h4># Maximum number of backtrack steps.</h4></html>", 0, 10, 1, Orientation.HORIZONTAL);
		backtrackMax.setPreferredSize(new Dimension(700, 20));
		backtrackMax.setMaximumSize(new Dimension(700, 20));
		backtrackMax.setMinimumSize(new Dimension(700, 20));
		add(backtrackMax, "0, 2, c, t");

		// max backtrackThreshold
		backtrackThreshold = slickerFactoryInstance.createNiceDoubleSlider(
				"<html><h4># Backtracking threshold.</h4></html>", 0, 4.0, 2.0, Orientation.HORIZONTAL);
		backtrackThreshold.setPreferredSize(new Dimension(700, 20));
		backtrackThreshold.setMaximumSize(new Dimension(700, 20));
		backtrackThreshold.setMinimumSize(new Dimension(700, 20));
		add(backtrackThreshold, "0, 3, c, t");

		//TODO: Add choice of Solver?

		setupUI(transCol, evClassCol, defMoveModelCost, defSyncCost, defMoveLogCost, slickerFactoryInstance, 4);

	}

	public int getExpectedModelMoves() {
		return modelMoves.getValue();
	}

	public int getBacktrackLimit() {
		return backtrackMax.getValue();
	}

	public double getBacktrackThreshold() {
		return backtrackThreshold.getValue();
	}

}
