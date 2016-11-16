package org.processmining.antialignments.test;

import java.util.Collection;
import java.util.Random;

import org.deckfour.xes.extension.std.XConceptExtension;
import org.deckfour.xes.factory.XFactoryRegistry;
import org.deckfour.xes.model.XEvent;
import org.deckfour.xes.model.XLog;
import org.deckfour.xes.model.XTrace;
import org.processmining.contexts.uitopia.UIPluginContext;
import org.processmining.contexts.uitopia.annotations.UITopiaVariant;
import org.processmining.framework.connections.ConnectionCannotBeObtained;
import org.processmining.framework.plugin.annotations.Plugin;
import org.processmining.framework.plugin.annotations.PluginLevel;
import org.processmining.framework.plugin.annotations.PluginVariant;
import org.processmining.models.connections.petrinets.behavioral.InitialMarkingConnection;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

@Plugin(name = "Simulate Petrinet with noise", level = PluginLevel.Local, //
returnLabels = { "log" }, returnTypes = { XLog.class },//
parameterLabels = { "Petrinet" }, //
help = "Simulate Petrinet with increasing noise level.", userAccessible = true)
public class PetriNetSimulator {

	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	pack = "AntiAlignments")
	@PluginVariant(variantLabel = "default", requiredParameterLabels = { 0 })
	public XLog export(UIPluginContext context, Petrinet net) throws ConnectionCannotBeObtained,
			IllegalTransitionException {

		XLog log = XFactoryRegistry.instance().currentDefault().createLog();

		InitialMarkingConnection initCon = context.getConnectionManager().getFirstConnection(
				InitialMarkingConnection.class, context, net);
		Marking initialMarking = initCon.getObjectWithRole(InitialMarkingConnection.MARKING);

		Random random = new Random(542323424643l);
		Transition[] trans = (Transition[]) net.getTransitions().toArray();

		//generate 25 traces
		for (int i = 0; i < 25; i++) {
			XTrace trace = XFactoryRegistry.instance().currentDefault().createTrace();
			log.add(trace);

			PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(net.getClass());
			semantics.setCurrentState(initialMarking);
			Collection<Transition> enabled = semantics.getExecutableTransitions();
			while (!enabled.isEmpty()) {
				Transition t = enabled.iterator().next();
				semantics.executeExecutableTransition(t);

				if (random.nextInt(25) > i) {
					// no noise
					XEvent event = XFactoryRegistry.instance().currentDefault().createEvent();
					XConceptExtension.instance().assignName(event, t.getLabel());
					trace.add(event);
				} else {
					// noise
					double d = random.nextDouble();
					if (d < 0.33) {
						// random name for event
						XEvent event = XFactoryRegistry.instance().currentDefault().createEvent();
						t = trans[random.nextInt(trans.length)];
						XConceptExtension.instance().assignName(event, t.getLabel());
						trace.add(event);
					} else if (d < 0.66) {
						// insert add extra event
						XEvent event = XFactoryRegistry.instance().currentDefault().createEvent();
						t = trans[random.nextInt(trans.length)];
						XConceptExtension.instance().assignName(event, t.getLabel());
						trace.add(event);

						event = XFactoryRegistry.instance().currentDefault().createEvent();
						t = trans[random.nextInt(trans.length)];
						XConceptExtension.instance().assignName(event, t.getLabel());
						trace.add(event);
					} else {
						// ignore event
					}
				}

				enabled = semantics.getExecutableTransitions();
			}

		}

		return log;
	}
}
