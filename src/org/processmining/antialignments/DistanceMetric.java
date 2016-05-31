package org.processmining.antialignments;

public interface DistanceMetric {

	public int getDistance(short[] trace, short[] word);

	public static class Hamming implements DistanceMetric {

		// simple distance based on indices
		public int getDistance(short[] trace, short[] word) {
			int d = 0;
			int i;
			for (i = 0; i < trace.length && i < word.length; i++) {
				if (trace[i] != word[i]) {
					d++;
				}
			}
			return d;
		}
	}

	public static class Edit implements DistanceMetric {

		// true edit distance
		public int getDistance(short[] trace, short[] word) {
			int len1 = trace.length;
			int len2 = word.length;

			// len1+1, len2+1, because finally return dp[len1][len2]
			int[][] dp = new int[len1 + 1][len2 + 1];

			for (int i = 0; i <= len1; i++) {
				dp[i][0] = i;
			}

			for (int j = 0; j <= len2; j++) {
				dp[0][j] = j;
			}

			// iterate though, and check last char
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {

					// if last two chars equal
					if (trace[i] == word[j]) {
						// update dp value for +1 length
						dp[i + 1][j + 1] = dp[i][j];
					} else {

						int replace = dp[i][j] + 1;
						int insert = dp[i][j + 1] + 1;
						int delete = dp[i + 1][j] + 1;

						int min = replace > insert ? insert : replace;
						min = delete > min ? min : delete;
						dp[i + 1][j + 1] = min;
					}
				}
			}

			return dp[len1][len2];
		}
	}

	public static class LongestCommonSubsequence implements DistanceMetric {

		// true edit distance
		public int getDistance(short[] trace, short[] word) {
			int len1 = trace.length;
			int len2 = word.length;

			// len1+1, len2+1, because finally return dp[len1][len2]
			int[][] dp = new int[len1 + 1][len2 + 1];

			for (int i = 0; i <= len1; i++) {
				dp[i][0] = i;
			}

			for (int j = 0; j <= len2; j++) {
				dp[0][j] = j;
			}

			// iterate though, and check last char
			for (int i = 0; i < len1; i++) {
				for (int j = 0; j < len2; j++) {

					// if last two chars equal
					if (trace[i] == word[j]) {
						// update dp value for +1 length
						dp[i + 1][j + 1] = dp[i][j];
					} else {

						// int replace = dp[i][j] + 1;
						int insert = dp[i][j + 1] + 1;
						int delete = dp[i + 1][j] + 1;

						dp[i + 1][j + 1] = insert > delete ? delete : insert;
					}
				}
			}

			return dp[len1][len2];
		}
	}

}
