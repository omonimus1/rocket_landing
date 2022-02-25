package coursework;

import model.Fitness;
import model.Individual;
import model.NeuralNetwork;

public class ExampleHillClimber extends NeuralNetwork{

	public static void main(String[] args) {
		NeuralNetwork hillClimber = new ExampleHillClimber();
		hillClimber.run();
	}
	
	@Override
	public void run() {
		//initialise a single individual
		best = new Individual();
		
		//run for max evaluations
		for(int gen = 0; gen < Parameters.maxEvaluations; gen++) {
			//mutate the best
			Individual candidate = mutateBest();
			
			//accept if better
			if(candidate.fitness < best.fitness) {
				best = candidate;
			}
			
			outputStats();
		}
		saveNeuralNetwork();
	}

	private Individual mutateBest() {
		Individual candidate = best.copy();
		for (int i = 0; i < candidate.chromosome.length; i++) {
			if (Parameters.random.nextDouble() < Parameters.mutateRate) {
				if (Parameters.random.nextBoolean()) {
					candidate.chromosome[i] += (Parameters.mutateChange);
				} else {
					candidate.chromosome[i] -= (Parameters.mutateChange);
				}
			}
		}
		Fitness.evaluate(candidate, this);
		return candidate;
	}

	
	@Override
	public double activationFunction(double x) {
		if (x < -20.0) {
			return -1.0;
		} else if (x > 20.0) {
			return 1.0;
		}
		return Math.tanh(x);
	}

	
}
