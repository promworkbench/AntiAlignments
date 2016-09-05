package org.processmining.antialignments;

import nl.tue.astar.AStarException;

import org.deckfour.xes.extension.std.XConceptExtension;
import org.deckfour.xes.model.XLog;
import org.processmining.antialignments.ilp.antialignment.AntiAlignmentValues;
import org.processmining.antialignments.ilp.antialignment.HeuristicAntiAlignmentAlgorithm;
import org.processmining.antialignments.ilp.util.AntiAlignments;
import org.processmining.contexts.uitopia.UIPluginContext;
import org.processmining.contexts.uitopia.annotations.UITopiaVariant;
import org.processmining.framework.connections.ConnectionCannotBeObtained;
import org.processmining.framework.connections.ConnectionManager;
import org.processmining.framework.plugin.PluginContext;
import org.processmining.framework.plugin.annotations.Plugin;
import org.processmining.framework.plugin.annotations.PluginLevel;
import org.processmining.framework.plugin.annotations.PluginVariant;
import org.processmining.models.connections.petrinets.EvClassLogPetrinetConnection;
import org.processmining.models.connections.petrinets.PNRepResultAllRequiredParamConnection;
import org.processmining.models.connections.petrinets.behavioral.FinalMarkingConnection;
import org.processmining.models.connections.petrinets.behavioral.InitialMarkingConnection;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;
import org.processmining.plugins.petrinet.replayer.PNLogReplayer;
import org.processmining.plugins.petrinet.replayresult.PNRepResult;

@Plugin(name = "Anti-Alignment Precision/Generalization", level = PluginLevel.NightlyBuild, //
returnLabels = { "Anti-alignments" }, returnTypes = { PNRepResult.class },//
parameterLabels = { "Petri net", "Event Log", "Alignment" }, //
help = "Measure precision/generalization using anti alignments.", userAccessible = true)
public class AntiAlignmentPlugin {

	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	pack = "AntiAlignments")
	@PluginVariant(variantLabel = "Without alignments", requiredParameterLabels = { 0, 1 })
	public PNRepResult measurePrecisionWithoutAlignments(final UIPluginContext context, Petrinet net, XLog log) {
		// replay log on model (or obtain existing replay result)
		PNRepResult alignments;
		try {
			PNLogReplayer replayer = new PNLogReplayer();
			alignments = replayer.replayLog(context, net, log);

			context.getProvidedObjectManager().createProvidedObject("Petrinet log replay result", alignments, context);

			return measurePrecision(context, net, log, alignments);

		} catch (ConnectionCannotBeObtained e) {
			context.log("No connection between the given net and log. For computing anti-alignment based precision, the plugin"
					+ " needs a Petri net with an initial and final marking!");
		} catch (AStarException e) {
			context.log("Error in replay algorithm: " + e.getMessage());
		}

		return null;
	}

	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	pack = "AntiAlignments")
	@PluginVariant(variantLabel = "With alignments", requiredParameterLabels = { 0, 1, 2 })
	public PNRepResult measurePrecision(PluginContext context, Petrinet net, XLog log, PNRepResult alignments) {
		// retrieve mapping between process model to log
		try {
			ConnectionManager connManager = context.getConnectionManager();
			EvClassLogPetrinetConnection conn = connManager.getFirstConnection(EvClassLogPetrinetConnection.class,
					context, net, log);
			TransEvClassMapping mapping = (TransEvClassMapping) conn
					.getObjectWithRole(EvClassLogPetrinetConnection.TRANS2EVCLASSMAPPING);

			// get marking
			InitialMarkingConnection initMarkingConn = connManager.getFirstConnection(InitialMarkingConnection.class,
					context, net);
			Marking initialMarking = initMarkingConn.getObjectWithRole(InitialMarkingConnection.MARKING);

			FinalMarkingConnection finalMarkingConn = connManager.getFirstConnection(FinalMarkingConnection.class,
					context, net);
			Marking finalMarking = finalMarkingConn.getObjectWithRole(FinalMarkingConnection.MARKING);

			PNRepResult replayRes = basicCodeStructureWithAlignments(net, initialMarking, finalMarking, log,
					alignments, mapping);
			if (replayRes != null) {
				context.addConnection(new PNRepResultAllRequiredParamConnection("Connection between replay result, "
						+ XConceptExtension.instance().extractName(log) + ", and " + net.getLabel(), net, log, mapping,
						null, null, replayRes));

			}
			return replayRes;

		} catch (ConnectionCannotBeObtained noConnection) {
			context.log("No connection between the given net and log. For computing anti-alignment based precision, the plugin"
					+ " needs a Petri net with an initial and final marking!");
		}

		return null;
	}

	public PNRepResult basicCodeStructureWithAlignments(Petrinet net, Marking initialMarking, Marking finalMarking,
			XLog xLog, PNRepResult alignments, TransEvClassMapping mapping) {

		HeuristicAntiAlignmentAlgorithm algorithm = new HeuristicAntiAlignmentAlgorithm(net, initialMarking,
				finalMarking, xLog, alignments, mapping);

		AntiAlignments aa = algorithm.computeAntiAlignments();

		AntiAlignmentValues values = algorithm.computePrecisionAndGeneralization(aa);

		return algorithm.getPNRepResult(aa, values);

	}

}
