package org.processmining.antialignments;

import java.text.NumberFormat;

import org.processmining.framework.util.HTMLToString;

public class AntiAlignmentValues implements HTMLToString {

	private final double alpha;
	private final double fitness;
	private final double precisionLog;
	private final double precisionTrace;
	private final double generalizationLog;
	private final double generalizationTrace;

	public AntiAlignmentValues(double alpha, double fitness, double precisionTrace, double precisionLog,
			double generalizationTrace, double generalizationLog) {
		this.alpha = alpha;
		this.fitness = fitness;
		this.precisionLog = precisionLog;
		this.precisionTrace = precisionTrace;
		this.generalizationLog = generalizationLog;
		this.generalizationTrace = generalizationTrace;

	}

	public double getFitness() {
		return fitness;
	}

	public double getPrecision() {
		return alpha * precisionTrace + (1 - alpha) * precisionLog;
	}

	public double getGeneralization() {
		return alpha * generalizationTrace + (1 - alpha) * generalizationLog;
	}

	public String toHTMLString(boolean includeHTMLTags) {
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(5);
		nf.setMaximumFractionDigits(5);

		StringBuffer buffer = new StringBuffer();

		if (includeHTMLTags) {
			buffer.append("<html>");
		}
		buffer.append("<head>");
		buffer.append("<title>Anti-Alignment-based Conformance checking</title>");
		buffer.append("</head><body><font fact=\"Arial\">");

		buffer.append("<h2>Fitness : " + nf.format(getFitness()) + "</h2>");
		buffer.append("<h2>Precision : " + nf.format(getPrecision()) + "</h2>");
		buffer.append("<h4>   Trace-Precision : " + nf.format(precisionTrace) + "</h4>");
		buffer.append("<h4>   Log-Precision  : " + nf.format(precisionLog) + "</h4>");
		buffer.append("<h2>Generalization : " + nf.format(getGeneralization()) + "</h2>");
		buffer.append("<h4>   Trace-Generalization : " + nf.format(generalizationTrace) + "</h4>");
		buffer.append("<h4>   Log-Generalization  : " + nf.format(generalizationLog) + "</h4>");

		buffer.append("</font></body>");
		if (includeHTMLTags) {
			buffer.append("</html>");
		}
		return buffer.toString();
	}

	public double getTracePrecision() {
		return precisionTrace;
	}

	public double getLogPrecision() {
		return precisionLog;
	}

	public double getTraceGeneralization() {
		return generalizationTrace;
	}

	public double getLogGeneralization() {
		return generalizationLog;
	}
}
