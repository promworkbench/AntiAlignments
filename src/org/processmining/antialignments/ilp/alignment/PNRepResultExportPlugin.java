package org.processmining.antialignments.ilp.alignment;

import java.io.PrintStream;

import org.processmining.plugins.petrinet.replayresult.PNRepResult;
import org.processmining.plugins.replayer.replayresult.SyncReplayResult;

// @Plugin(name = "Write PN rep result", level = PluginLevel.Local, //
// returnLabels = { "nothing" }, returnTypes = { Object.class },//
// parameterLabels = { "Replay Result" }, //
// help = "Measure precision/generalization using anti alignments.",
// userAccessible = true)
public class PNRepResultExportPlugin {

	//	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	//	pack = "AntiAlignments")
	//	@PluginVariant(variantLabel = "default", requiredParameterLabels = { 0 })
	//	public Object export(UIPluginContext context, PNRepResult result) {
	//		PrintStream out;
	//		try {
	//			out = new PrintStream(new File("d:/temp/antialignment/exported.csv"));
	//		} catch (FileNotFoundException e1) {
	//			out = System.out;
	//		}
	//		String sep = ";";
	//
	//		printResult(out, null, 0, 0, sep);
	//
	//		printResult(out, result, -1, -1, sep);
	//
	//		out.close();
	//		return null;
	//
	//	}

	public void printResult(PrintStream out, PNRepResult result, int estRows, int estColumns, String sep) {
		if (result == null) {
			out.print("Solver");
			out.print(sep);
			out.print("FirstTrace");
			out.print(sep);
			out.print("RepresentedCount");
			out.print(sep);
			out.print("XLength");
			out.print(sep);
			out.print("minEvent");
			out.print(sep);
			out.print("BacktrackBound");
			out.print(sep);
			out.print("BacktrackThreshold");
			out.print(sep);
			out.print("FitnessCost");
			out.print(sep);
			out.print("TotalFitnessCost");
			out.print(sep);
			out.print("MoveLogFitness");
			out.print(sep);
			out.print("MoveModelFitness");
			out.print(sep);
			out.print("SolveSteps");
			out.print(sep);
			out.print("TraceFitness");
			out.print(sep);
			out.print("Time");
			out.print(sep);
			out.print("OriginalTraceLength");
			out.print(sep);
			out.print("AlignmentLength");
			out.print(sep);
			out.print("Reliable");
			out.print(sep);
			out.print("EstimatedRows");
			out.print(sep);
			out.print("EstimatedColumns");
			out.println();

		} else {
			for (SyncReplayResult r : result) {
				out.print(result.getInfo().get(HeuristicPNetReplayerAlgorithm.SOLVER));
				out.print(sep);
				out.print(r.getTraceIndex().first());
				out.print(sep);
				out.print(r.getTraceIndex().size());
				out.print(sep);
				out.print(result.getInfo().get(HeuristicPNetReplayerAlgorithm.CUTOFF));
				out.print(sep);
				out.print(result.getInfo().get(HeuristicPNetReplayerAlgorithm.MINEVENT));
				out.print(sep);
				out.print(result.getInfo().get(HeuristicPNetReplayerAlgorithm.BACKTRACKBOUND));
				out.print(sep);
				out.print(result.getInfo().get(HeuristicPNetReplayerAlgorithm.BACKTRACKTHRESHOLD));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.RAWFITNESSCOST));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.RAWFITNESSCOST) * r.getTraceIndex().size());
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.MOVELOGFITNESS));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.MOVEMODELFITNESS));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.NUMSTATEGENERATED));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.TRACEFITNESS));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.TIME));
				out.print(sep);
				out.print(r.getInfo().get(PNRepResult.ORIGTRACELENGTH));
				out.print(sep);
				out.print(r.getStepTypes().size());
				out.print(sep);
				out.print(estRows);
				out.print(sep);
				out.print(estColumns);
				out.print(sep);
				out.print(Boolean.toString(r.isReliable()));
				out.println();
			}
		}
		out.flush();
	}
}
