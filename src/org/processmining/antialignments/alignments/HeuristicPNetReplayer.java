package org.processmining.antialignments.alignments;

import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gnu.trove.map.hash.TShortObjectHashMap;

import java.util.Collection;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.classification.XEventClassifier;
import org.deckfour.xes.model.XEvent;
import org.deckfour.xes.model.XLog;
import org.deckfour.xes.model.XTrace;
import org.processmining.antialignments.AntiAlignmentValues;
import org.processmining.contexts.uitopia.UIPluginContext;
import org.processmining.contexts.uitopia.annotations.UITopiaVariant;
import org.processmining.framework.connections.Connection;
import org.processmining.framework.connections.ConnectionCannotBeObtained;
import org.processmining.framework.connections.annotations.ConnectionObjectFactory;
import org.processmining.framework.plugin.PluginContext;
import org.processmining.framework.plugin.PluginExecutionResult;
import org.processmining.framework.plugin.PluginParameterBinding;
import org.processmining.framework.plugin.annotations.Plugin;
import org.processmining.framework.plugin.annotations.PluginLevel;
import org.processmining.framework.plugin.annotations.PluginVariant;
import org.processmining.framework.util.Pair;
import org.processmining.models.connections.petrinets.EvClassLogPetrinetConnection;
import org.processmining.models.connections.petrinets.behavioral.FinalMarkingConnection;
import org.processmining.models.connections.petrinets.behavioral.InitialMarkingConnection;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetGraph;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;

@Plugin(name = "Heuristic alignment calculator", level = PluginLevel.NightlyBuild, //
returnLabels = { "Approximate alignments" }, returnTypes = { AntiAlignmentValues.class },//
parameterLabels = { "Petri net", "Event Log" }, //
help = "Compute approximate alignments using fast heuristics.", userAccessible = true)
public class HeuristicPNetReplayer {

	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	pack = "AntiAlignments")
	@PluginVariant(variantLabel = "Heuristic", requiredParameterLabels = { 0, 1 })
	public void computeAlignments(UIPluginContext context, Petrinet net, XLog xLog) throws Exception {

		EvClassLogPetrinetConnection conn = null;
		Marking initialMarking = null;
		Marking finalMarking = null;

		// check existence of initial marking
		try {
			InitialMarkingConnection initCon = context.getConnectionManager().getFirstConnection(
					InitialMarkingConnection.class, context, net);
			initialMarking = (Marking) initCon.getObjectWithRole(InitialMarkingConnection.MARKING);
			if (initialMarking.isEmpty()) {
				JOptionPane
						.showMessageDialog(
								new JPanel(),
								"The initial marking is an empty marking. If this is not intended, remove the currently existing InitialMarkingConnection object and then use \"Create Initial Marking\" plugin to create a non-empty initial marking.",
								"Empty Initial Marking", JOptionPane.INFORMATION_MESSAGE);
			}
		} catch (ConnectionCannotBeObtained exc) {
			if (JOptionPane.YES_OPTION == JOptionPane.showConfirmDialog(new JPanel(),
					"No initial marking is found for this model. Do you want to create one?", "No Initial Marking",
					JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE)) {
				initialMarking = createMarking(context, net, InitialMarkingConnection.class);
			}
		}

		// check existence of final marking
		try {
			FinalMarkingConnection finalCon = context.getConnectionManager().getFirstConnection(
					FinalMarkingConnection.class, context, net);
			finalMarking = finalCon.getObjectWithRole(FinalMarkingConnection.MARKING);
		} catch (ConnectionCannotBeObtained exc) {
			if (JOptionPane.YES_OPTION == JOptionPane.showConfirmDialog(new JPanel(),
					"No final marking is found for this model. Do you want to create one?", "No Final Marking",
					JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE)) {
				finalMarking = createMarking(context, net, FinalMarkingConnection.class);
			}
		}

		// check connection in order to determine whether mapping step is needed
		// of not
		try {
			// connection is found, no need for mapping step
			// connection is not found, another plugin to create such connection
			// is automatically
			// executed
			conn = context.getConnectionManager().getFirstConnection(EvClassLogPetrinetConnection.class, context, net,
					xLog);
		} catch (Exception e) {
			JOptionPane.showMessageDialog(new JPanel(), "No mapping can be constructed between the net and the log");
			return;
		}

		// So, we have setup everything
		// Now create the replayer.
		TObjectIntMap<TShortList> tempLog = new TObjectIntHashMap<>(xLog.size());
		TShortObjectMap<String> short2label = new TShortObjectHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);
		TObjectShortMap<String> label2short = new TObjectShortHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);
		
		short c = 0;

		XEventClassifier classifier = (XEventClassifier) conn
				.getObjectWithRole(EvClassLogPetrinetConnection.EVENTCLASSIFIER);
		TransEvClassMapping mapping = (TransEvClassMapping) conn
				.getObjectWithRole(EvClassLogPetrinetConnection.TRANS2EVCLASSMAPPING);

		for (XTrace trace : xLog) {
			TShortList list = new TShortArrayList(trace.size());
			for (XEvent event : trace) {
				String clazz = classifier.getClassIdentity(event);
				short id = label2short.putIfAbsent(clazz, c);
				if (id == label2short.getNoEntryValue()) {
					short2label.put(c, clazz);
					id = c;
					c++;
				}
				list.add(id);
			}
			tempLog.adjustOrPutValue(list, 1, 1);
		}
		for (Transition t : net.getTransitions()) {
			XEventClass clazz = mapping.get(t);
			short id = label2short.putIfAbsent(clazz.getId(), c);
			if (id == label2short.getNoEntryValue()) {
				short2label.put(c, clazz.getId());
				id = c;
				c++;
			}
			label2short.putIfAbsent(t.getLabel(), id);
		}

		int t = 0;
		int[] frequencies = new int[tempLog.size()];
		short[][] log = new short[tempLog.size()][];
		TObjectIntIterator<TShortList> it = tempLog.iterator();
		while (it.hasNext()) {
			it.advance();
			log[t] = it.key().toArray();
			frequencies[t] = it.value();
			t++;
		}

		AlignmentILPCalculator calculator = new AlignmentILPCalculator(net, initialMarking, finalMarking, label2short,
				short2label, log);
		calculator.solveSequential(initialMarking, finalMarking, 0);

	}

	private Marking createMarking(UIPluginContext context, PetrinetGraph net, Class<? extends Connection> classType) {
		Collection<Pair<Integer, PluginParameterBinding>> plugins = context.getPluginManager().find(
				ConnectionObjectFactory.class, classType, context.getClass(), true, false, false, net.getClass());
		PluginContext c2 = context.createChildContext("Creating connection of Type " + classType);
		Pair<Integer, PluginParameterBinding> pair = plugins.iterator().next();
		PluginParameterBinding binding = pair.getSecond();
		try {
			PluginExecutionResult pluginResult = binding.invoke(c2, net);
			pluginResult.synchronize();
			context.getProvidedObjectManager().createProvidedObjects(c2); // push the objects to main context
			return (Marking) pluginResult.getResults()[pair.getFirst()];
		} catch (Exception e) {
		} finally {
			c2.getParentContext().deleteChild(c2);
		}
		return null;
	}

}
