import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

public class picomap {

	public static void main(String[] args) {
		ObjectInputStream in = null;
		String genome_in = null;
		FM fm = null;
		BufferedReader rin = null;
		ArrayList<String[]> readsIn = new ArrayList<String[]>();
		String[][] reads;
		String[] read = new String[2];
		String[] finalTotalResult;
		int resultIndex = 0;
		
		try {
			in = new ObjectInputStream(new FileInputStream((args[0])));
			fm = (FM) in.readObject();
			
			byte[] genomeBytes = new byte[fm.SA.length];
			in.readFully(genomeBytes);
			genome_in = new String(genomeBytes);
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			try {
				if (in != null) {
					in.close();
				}
			} catch (IOException e) {
				System.err.println("Couldn't close " + args[0]);
			}
		}
		
		try {
			rin = new BufferedReader(new FileReader(args[1]));
			String line; 
			while((line = rin.readLine()) != null) {
				read[0] = line;
				read[1] = rin.readLine();
				readsIn.add(read.clone());
			}
			
			reads = readsIn.toArray(new String[0][0]);
			finalTotalResult = new String[reads.length];
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			try {
				if (in != null) {
					rin.close();
				}
			} catch (IOException e) {
				System.err.println("Couldn't close " + args[1]);
			}
		}
		
		int mismatch = Integer.valueOf(args[2]);
		int gap = Integer.valueOf(args[3]);
		int[] results;
		
		for (int i = 0; i < reads.length; i++) {
			results = completeQuery(fm, reads[i][1]);
			finalTotalResult[resultIndex] = reads[i][0].substring(1) + "\t" + results.length;
			for (int j = 0; j < results.length; j++) {
				finalTotalResult[resultIndex] += "\n" + results[j] + "\t0\t" + reads[i][1].length() + "=";
			}
			if (results.length == 0) {
				//project 5 all reads are of length 100
				//break into 5 reads starting at intervals of 20 (length 100 / 5) and all ending at 100
				
				int section = reads[i][1].length() / 5;
				int curSeed = 0;
				int score = 0;
				int[][] subSeeds = new int[5][];
				ArrayList<int[]> scores = new ArrayList<int[]>();
				ArrayList<int[]> best = new ArrayList<int[]>();
				boolean found = false;
				
				for (int j = 0; j < 5; j++) subSeeds[j] = partialQuery(fm, 
						reads[i][1].substring(section * j, section * (j + 1)));
				for (int seedx = 0; seedx < 5; seedx++) {
					for (int seedy = 0; seedy < subSeeds[seedx].length; seedy++) {
						curSeed = subSeeds[seedx][seedy];
						score = 0;
						for (int ssx = seedx + 1; ssx < 5; ssx++) {
							found = false;
							for (int ssy = 0; ssy < subSeeds[ssx].length && !found; ssy++) {
								if (curSeed + (section * (ssx - seedx)) - 15 < subSeeds[ssx][ssy] &&
										curSeed + (section * (ssx - seedx)) + 15 > subSeeds[ssx][ssy]) {
									score++;
								} else if (curSeed + (section * (ssx - seedx)) + 15 <= subSeeds[ssx][ssy]) found = true;
							}
						}
						found = false;
						if (score > 0) {
						int[] finalScore = new int[3];
							finalScore[0] = subSeeds[seedx][seedy];
							finalScore[1] = score;
							finalScore[2] = seedx;
							scores.add(finalScore);
						}
					}
				}
				score = 0;
				for (int[] curScore : scores) if (curScore[1] >= score) score = curScore[1];
				for (int[] curScore : scores) {
					if (curScore[1] > score) System.out.println("scoring error");
					else if (curScore[1] == score) best.add(curScore);
					for (int j = curScore[2]; j > 0; j--) curScore[0] -= section;
					curScore[0] -= 15;
				}
				Object[][] alignResults = new Object[best.size()][3];
				score = Integer.MAX_VALUE;
				found = true;
				
				//fittingAlign: start, score, cigar
				for (int j = 0; j < best.size(); j++) {
					int start = best.get(j)[0];
					alignResults[j] = fittingAlign(reads[i][1],
							genome_in.substring(start, start + reads[i][1].length() + 30), mismatch, gap);
					alignResults[j][0] = ((Integer) alignResults[j][0]) + start;
					if ((Integer) alignResults[j][1] <= score) {
						score = (Integer) alignResults[j][1];
					}
				}
				HashMap<Integer, Object[]> seeds = new HashMap<Integer, Object[]>();
				ArrayList<Integer> mine = new ArrayList<Integer>();
				for (int j = 0; j < alignResults.length; j++) {
					if ((Integer) alignResults[j][1] < score) System.out.println("scoring issue");
					else if ((Integer) alignResults[j][1] == score) {
						if (seeds.get(alignResults[j][0]) == null) {
							seeds.put((Integer) alignResults[j][0], alignResults[j]);
							mine.add((Integer) alignResults[j][0]);
						}
					} 
				}
				//Collections.sort(mine, Collections.reverseOrder());
				Collections.sort(mine);
				for (Integer current : mine) {
					Object[] currentArray = seeds.get(current);
					if (found) {
						finalTotalResult[resultIndex] = reads[i][0].substring(1) + "\t" + mine.size();
						found = false;
					}
					finalTotalResult[resultIndex] += "\n" + ((Integer) currentArray[0]) + "\t-" + (Integer) currentArray[1] + "\t" +  currentArray[2];
				}
			}
			resultIndex++;
		}
		
		try {
			FileWriter outStream = new FileWriter(new File(args[4]));
			
			for (int i = 0; i < finalTotalResult.length; i++) {
				if (i != 0) outStream.write("\n");
				outStream.write(finalTotalResult[i]);
			}
			
			outStream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static int[] completeQuery(FM fm, String query) {
		char[] nucs;
		int[] startingPoses = new int[5];
		startingPoses[0] = 0;
		for (int i = 1; i < 5; i++) startingPoses[i] = startingPoses[i - 1] + fm.first[i - 1];
		int letter, matchLen, startIndex, endIndex, numHits;
		
		
		numHits = Integer.MAX_VALUE;
		nucs = query.toCharArray();
		matchLen = 1;
		letter = parNum(nucs[nucs.length - 1]);
		startIndex = startingPoses[letter];
		endIndex = startIndex + fm.first[letter] - 1;
		for (int j = nucs.length - 2; numHits > 0 && j > -1; j--) {
			letter = parNum(nucs[j]);
			numHits = fm.tally[letter][endIndex] - fm.tally[letter][startIndex - 1];
			startIndex = startingPoses[letter] + fm.tally[letter][startIndex - 1];
			endIndex = startIndex + numHits - 1;
			if (numHits > 0) matchLen++; 
		}
		if (matchLen < nucs.length) return new int[0];
		else {
			int[] result = new int[numHits];
			for (int j = startIndex; j < startIndex + numHits; j++) {
				result[j - startIndex] = fm.SA[j];
			}
			Arrays.sort(result);
			return result;
		}
	}
	
	private static int[] partialQuery(FM fm, String query) {
		char[] nucs;
		int[] startingPoses = new int[5];
		startingPoses[0] = 0;
		for (int i = 1; i < 5; i++) startingPoses[i] = startingPoses[i - 1] + fm.first[i - 1];
		int letter, startIndex, endIndex, numHits, oldNumHits, oldIndex;
		numHits = Integer.MAX_VALUE;
		nucs = query.toCharArray();
		
		oldNumHits = oldIndex = 1;
		letter = parNum(nucs[nucs.length - 1]);
		startIndex = startingPoses[letter];
		endIndex = startIndex + fm.first[letter] - 1;
		for (int j = nucs.length - 2; numHits > 0 && j > -1; j--) {
			letter = parNum(nucs[j]);
			oldNumHits = numHits;
			oldIndex = startIndex;
			numHits = fm.tally[letter][endIndex] - fm.tally[letter][startIndex - 1];
			startIndex = startingPoses[letter] + fm.tally[letter][startIndex - 1];
			endIndex = startIndex + numHits - 1;
		}
		if (numHits == 0) {
			int[] result = new int[oldNumHits];
			for (int j = oldIndex; j < oldIndex + oldNumHits; j++) {
				result[j - oldIndex] = fm.SA[j];
			}
			Arrays.sort(result);
			return result;
		}
		else {
			int[] result = new int[numHits];
			for (int j = startIndex; j < startIndex + numHits; j++) {
				result[j - startIndex] = fm.SA[j];
			}
			Arrays.sort(result);
			return result;
		}
	}
	
	private static Object[] fittingAlign(String x, String y, int mismatch, int gap) {
		char[] query = x.toCharArray();
		char[] reference = y.toCharArray();
		Cell[][] vals = new Cell[reference.length + 1][query.length + 1];
		Direction[] path = new Direction[query.length + reference.length + 2];
		int pathLength = 0;
		
		for (int i = 1; i < vals.length; i++) { vals[i][0] = new Cell(0, Direction.LEFT); }
		for (int i = 1; i < vals[0].length; i++) { vals[0][i] = new Cell(gap * i, Direction.DOWN); }
		vals[0][0] = new Cell(0, Direction.NONE);
		
		for (int j = 1; j < vals[0].length; j++) {
			for (int i = 1; i < vals.length; i++) {
				if (reference[i - 1] == query[j - 1]) 
					vals[i][j] = new Cell(vals[i - 1][j - 1].getValue(), Direction.DIAG); 
				else {
					vals[i][j] = minVal(vals[i - 1][j], vals[i - 1][j - 1], vals[i][j - 1], mismatch, gap);
				}
			}
		}
		
		//fourth line: score, reference (y, second string) start, reference end,  CIGAR
		//CIGAR: 
		int integer1 = 0;
		int max = query.length;
		int start = 0;
		int end = 0;
		int maxVal = Integer.MAX_VALUE;
		int[] theseVals = new int[vals.length];
		for (int m = 0; m < vals.length; m++) {
			theseVals[m] = vals[m][max].getValue();
			if (vals[m][max].getValue() < maxVal) {
				maxVal = vals[m][max].getValue();
				integer1 = m;
			}
		}
		
		/*for (int j = vals[0].length - 1; j > -1; j--) {
			for (int i = 0; i < vals.length - 1; i++) {
				System.out.print(vals[i][j]);
			}
			System.out.println();
		}*/

		end = integer1;
		boolean hitEnd = false;
		for (Cell cur = vals[integer1][max]; cur != null && cur.getPointer() != null && max > -1; pathLength++) {
			if (max == 0 && !hitEnd) {
				hitEnd = true;
				start = integer1;
			} else {
				path[pathLength] = cur.getPointer();
				if (path[pathLength] == Direction.LEFT) cur = vals[--integer1][max];
				else if (path[pathLength] == Direction.DIAG) cur = vals[--integer1][--max];
				else if (path[pathLength] == Direction.MISD) cur = vals[--integer1][--max];
				else if (path[pathLength] == Direction.DOWN) cur = vals[integer1][--max];
				else if (path[pathLength] == Direction.NONE) cur = null;
			}
		}
		
		Object[] record = new Object[3];
		record[0] = start;
		record[1] = maxVal;
		Direction prev = null;
		integer1 = 0;
		String cigar = "";
		for (Direction curr : path) { if (curr != null) {
			if (curr == prev) integer1++;
			else {
				if (prev == Direction.LEFT) cigar = integer1 + "D" + cigar;
				else if (prev == Direction.DOWN) cigar = integer1 + "I" + cigar;
				else if (prev == Direction.DIAG) cigar = integer1 + "=" + cigar;
				else if (prev == Direction.MISD) cigar = integer1 + "X" + cigar;
				integer1 = 1;
			}
			
			prev = curr;
		}}
		if (start != 0) cigar = cigar.split("D", 2)[1];
		record[2] = cigar;
		return record;
	}
	
	enum Direction { DOWN, DIAG, MISD, LEFT, NONE}
	
	private static Cell minVal(Cell left, Cell diag, Cell down, int mis, int gap) {
		int leftVal = left.getValue() + gap;
		int diagVal = diag.getValue() + mis;
		int downVal = down.getValue() + gap;
		
		if (leftVal < diagVal) {
			if (leftVal < downVal) return new Cell(leftVal, Direction.LEFT);
			else return new Cell(downVal, Direction.DOWN);
		} else {
			if (diagVal < downVal) return new Cell(diagVal, Direction.MISD);
			else return new Cell(downVal, Direction.DOWN);
		}
	}
	
	private static class Cell {
		private int value;
		private Direction pointer;
		
		public Cell(int value, Direction pointer) {
			this.value = value;
			this.pointer = pointer;
		}
		
		public String toString() {
			return "   " + value + "|" + pointer + "   ";
		}
		
		public int getValue() { return value; }
		public Direction getPointer() {return pointer; }
	}
	
	private static int parNum(char x) {
		if (x == '$') return 0;
		if (x == 'A') return 1;
		if (x == 'C') return 2;
		if (x == 'G') return 3;
		if (x == 'T') return 4;
		else {
			System.out.println("invalid char found");
			return -1;
		}
	}
}
