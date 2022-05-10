import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;

public class picoindex {

	public static void main(String[] args) {
		String input = args[0];
		String output = args[1];

		StringBuffer genome_buf = new StringBuffer();

		try(BufferedReader br = new BufferedReader(new FileReader(input))) {
			for(String line; (line = br.readLine()) != null; ) {
				if (!line.startsWith(">")) {
					genome_buf.append(line.toUpperCase().strip());
				}
			}
			// line is not visible here.
		} catch (IOException e) {
			e.printStackTrace();
		}
		genome_buf.append("$");
		String genome = genome_buf.toString();
		int[] sa = SuffixArray.constructSuffixArray(genome);
		FM fm = SAtoFM(sa, genome);
		
		ObjectOutputStream out = null;
		try {
			out = new ObjectOutputStream(new FileOutputStream(output));
			out.writeObject(fm);
			out.writeBytes(genome);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			try {
				if (out != null) {
					out.close();
				}
			} catch (IOException e) {
				System.err.println("Couldn't close " + output);
			}
		}

		// Below is how the serialized file would be
		// read in.
		/*
		ObjectInputStream in = null;
		String genome_in = null;
		FM fm_in = null;
		try {
			in = new ObjectInputStream(new FileInputStream(output));
			fm_in = (FM) in.readObject();
			if (!fm_in.equals(fm)) {
				System.out.println("suffix arrays don't match");
			}
			byte[] genomeBytes = new byte[fm_in.SA.length];
			in.readFully(genomeBytes);
			genome_in = new String(genomeBytes);
			if (!genome_in.equals(genome)) {
				System.out.println("genomes don't match");
			}
			System.out.println("all matched!");
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			try {
				if (in != null) {
					in.close();
				}
			} catch (IOException e) {
				System.err.println("Couldn't close " + output);
			}
		}
		*/
		
	}

	
	public static FM SAtoFM(int[] SA, String genome) {
		FM ans = new FM();
		int[] firstCount = new int[5];
		int[] BWT = new int[SA.length];
		int[][] tallyTable = new int[5][SA.length];
 		int temp;
		
		for (int i = 0; i < SA.length; i++) {
			temp = SA[i];
			
			//BWT (or L vector) of genome
			temp--;
			if (temp == -1) temp = SA.length - 1;
			temp = parNum(genome.charAt(temp));
			BWT[i] = temp;
			
			//tally table (for L)
			if (i != 0) for (int j = 0; j < 5; j++) tallyTable[j][i] = tallyTable[j][i - 1];
			tallyTable[temp][i]++;
		}
		
		//first count
		for (int i = 0; i < 5; i++) firstCount[i] = tallyTable[i][SA.length - 1];
		
		ans.SA = SA;
		ans.BWT = BWT;
		ans.first = firstCount;
		ans.tally = tallyTable;
		
		return ans;
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
