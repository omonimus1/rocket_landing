package coursework;

import java.lang.reflect.Field;
import java.util.Random;
import model.LunarParameters;
import model.NeuralNetwork;
import model.LunarParameters.DataSet;

public class Parameters {
	private static int numHidden = 12;
	private static int numGenes = calculateNumGenes();
	/* 
	 * minGene: specifies minimum and maximum weight values 
	 * Initial value: -1
	 */
	public static double minGene = -1; 
	public static double maxGene = +1; 
		
	public static int populationSize = 50; // number of list of solutions to crossover at each iteration; 
	public static int maxEvaluations = 20000; // cannot be set 20000
	
	//Random number generator used throughout the application
	public static long seed = System.currentTimeMillis();
	public static Random random = new Random(seed);

	//set the NeuralNetwork class here to use your code from the GUI
	public static Class neuralNetworkClass = ExampleEvolutionaryAlgorithm.class;
	
	public static int tournamentSize = 9; 
	
	public enum InitialisationType { POS_NEG, RANDOM,}
	public enum SelectionCategory { RANDOM_SELECTION, TOURNAMENT, ROULETTE_SELECTION, RANK_SELECTION}
	public enum CrossoverType { UNIFORM, ONE_POINT, TWO_POINTS}
	public enum MutationType { CONSTRAINED, STANDARD }
	public enum ReplaceType { TOURNAMENT, WORST  }
	public enum ActivationType { RELU, TANH, SELU}
	
	public static InitialisationType initialisationType = InitialisationType.POS_NEG;  
	public static SelectionCategory selectionCategory = SelectionCategory.TOURNAMENT;
	public static CrossoverType crossoverType = CrossoverType.TWO_POINTS;
	public static MutationType mutationType = MutationType.STANDARD; 
	public static ReplaceType replaceType = ReplaceType.TOURNAMENT;
	public static ActivationType activationType = ActivationType.SELU;
	
	
	// Rate = probability of changing a gene
	// Change = the maximum +/- adjustment to the gene value
	public static double mutateRate = 0.45; // Mutation rate for mutation operator
	public static double mutateChange = 0.95; // Delta change for mutation operator
	
	public static double SAcoolingRate = 0.0011;


	public static int getNumGenes() {					
		return numGenes;
	}
	
	private static int calculateNumGenes() {
		int num = (NeuralNetwork.numInput * numHidden) + (numHidden * NeuralNetwork.numOutput) + numHidden + NeuralNetwork.numOutput;
		return num;
	}

	public static int getNumHidden() {
		return numHidden;
	}
	
	public static void setHidden(int nHidden) {
		numHidden = nHidden;
		numGenes = calculateNumGenes();		
	}

	public static String printParams() {
		String str = "";
		for(Field field : Parameters.class.getDeclaredFields()) {
			String name = field.getName();
			Object val = null;
			try {
				val = field.get(null);
			} catch (Exception e) { 
				e.printStackTrace();
				throw new RuntimeException(e); 
			}
			str += name + " \t" + val + "\r\n";
			
		}
		return str;
	}
	
	public static void setDataSet(DataSet dataSet) {
		LunarParameters.setDataSet(dataSet);
	}
	
	public static DataSet getDataSet() {
		return LunarParameters.getDataSet();
	}
	
	public static void main(String[] args) {
		printParams();
	}
}
