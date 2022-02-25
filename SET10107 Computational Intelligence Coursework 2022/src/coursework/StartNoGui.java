package coursework;

import model.Fitness;
import model.LunarParameters.DataSet;
import model.NeuralNetwork;

/**
 * Example of how to to run the {@link ExampleEvolutionaryAlgorithm} without the need for the GUI
 * This allows you to conduct multiple runs programmatically 
 * The code runs faster when not required to update a user interface
 *
 */
public class StartNoGui {

	public static void main(String[] args) {
		/**
		 * Train the Neural Network using our Evolutionary Algorithm 
		 * 
		 */

		//Set the parameters here or directly in the Parameters Class
		Parameters.maxEvaluations = 20000; // Used to terminate the EA after this many generations
		Parameters.popSize = 200; // Population Size

		//number of hidden nodes in the neural network
		Parameters.setHidden(5);
		
		//Set the data set for training 
		Parameters.setDataSet(DataSet.Training);
		
		
		//Create a new Neural Network Trainer Using the above parameters 
		NeuralNetwork nn = new ExampleEvolutionaryAlgorithm();		
		
		//train the neural net (Go and have a coffee) 
		nn.run();
		
		/* Print out the best weights found
		 * (these will have been saved to disk in the project default directory using 
		 * the saveWeights method in EvolutionaryTrainer) 
		 */
		System.out.println(nn.best);
		
		
		
		
		/**
		 * The last File Saved to the Output Directory will contain the best weights /
		 * Parameters and Fitness on the Training Set 
		 * 
		 * We can used the trained NN to Test on the test Set
		 */
		Parameters.setDataSet(DataSet.Test);
		double fitness = Fitness.evaluate(nn);
		System.out.println("Fitness on " + Parameters.getDataSet() + " " + fitness);
		
		/**
		 * Or We can reload the NN from the file generated during training and test it on a data set 
		 * We can supply a filename or null to open a file dialog 
		 * Note that files must be in the project root and must be named *-n.txt
		 * where "n" is the number of hidden nodes
		 * ie  1518461386696-5.txt was saved at timestamp 1518461386696 and has 5 hidden nodes
		 * Files are saved automatically at the end of training
		 *  
		 */
		
		ExampleEvolutionaryAlgorithm nn2 = ExampleEvolutionaryAlgorithm.loadNeuralNetwork("1518446327913-5.txt");
		Parameters.setDataSet(DataSet.Random);
		double fitness2 = Fitness.evaluate(nn2);
		System.out.println("Fitness on " + Parameters.getDataSet() + " " + fitness2);
		
	}
}
