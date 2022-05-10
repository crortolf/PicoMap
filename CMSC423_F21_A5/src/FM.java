

public class FM implements java.io.Serializable{
	int[] SA, BWT, first;
	int[][] tally;
	
	public FM(int[] SA, int[] BWT, int[] first, int[][] tally) {
		this.SA = SA;
		this.BWT = BWT;
		this.first = first;
		this.tally = tally;
	}
	
	public FM() {
		SA = BWT = first = null;
		tally = null;
	}
	
	public boolean equals(FM other) {
		if (!java.util.Arrays.equals(SA, other.SA)) return false;
		if (!java.util.Arrays.equals(BWT, other.BWT)) return false;
		if (!java.util.Arrays.equals(first, other.first)) return false;
		if (!java.util.Arrays.deepEquals(tally, other.tally)) return false;
		return true;
	}
	}
