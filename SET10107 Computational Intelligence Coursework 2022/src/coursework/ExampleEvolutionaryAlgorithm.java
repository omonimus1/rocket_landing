package coursework;

import java.util.ArrayList;
import java.util.Collections;
import java.util.stream.Collector;
import java.util.stream.Collectors;

import coursework.Parameters.InitialisationType;
import coursework.Parameters.SelectionCategory;
import model.Fitness;
import model.Individual;
import model.NeuralNetwork;


public class ExampleEvolutionaryAlgorithm extends NeuralNetwork {
	// The Main Evolutionary Loop
	@Override
	public void run() {					
		// Initialise the population with pseudo-random weights
		if( Parameters.initialisationType == InitialisationType.POS_NEG)
			population = PosNegInitialise();
		else
			population = initialise();
		// Evolutionary Algorithm processing
		while (evaluations < Parameters.maxEvaluations) {
			
			// Select 2 Individuals from the current population
			Individual first_parent = null, second_parent = null;
			
			if( Parameters.selectionCategory == SelectionCategory.TOURNAMENT)
			{
				first_parent = tournamentSelect(); 
				second_parent = tournamentSelect();
			}
			else if( Parameters.selectionCategory == SelectionCategory.RANK_SELECTION)
			{
				first_parent = rankSelect(); 
				second_parent = rankSelect();
			}
			else if( Parameters.selectionCategory == SelectionCategory.ROULETTE_SELECTION)
			{
				first_parent = rouletteSelect(); 
				second_parent = rouletteSelect();
			}
			else if( Parameters.selectionCategory == SelectionCategory.RANDOM_SELECTION)
			{
				first_parent = randSelect(); 
				second_parent = randSelect();
			}
	
			// Generate children_set by crossover
			ArrayList<Individual> children_dataset;
			switch(Parameters.crossoverType) {
				case ONE_POINT:
					children_dataset = onePointCrossover(first_parent, second_parent);
					break;
				case UNIFORM:
					children_dataset = uniformCrossover(first_parent, second_parent);
					break;
				default:
					children_dataset = twoPointCrossover(first_parent, second_parent);
					break;
			}
				
			
			switch(Parameters.mutationType) {
				case CONSTRAINED:
					constrainedMutation(children_dataset);
					break;
				default:
					mutate(children_dataset);
					break;
			}
			
			evaluateIndividuals(children_dataset);
			// Replace children_set in population
			switch(Parameters.replaceType) {
				case TOURNAMENT:
				default:
					tournamentReplace(children_dataset);
					break;
				case WORST:
					replaceWorst(children_dataset);
					break;
			}
			
			best = getBest();
			// System.out.println("Best Individual obtained from the population: " + best);
			outputStats();
		}
		//save the trained network to disk
		saveNeuralNetwork();
	}


	/**
	 * Sets the fitness of the individuals passed as parameters (whole population)
	 */
	private void evaluateIndividuals(ArrayList<Individual> individuals) {
		for (Individual individual : individuals) {
			individual.fitness = Fitness.evaluate(individual, this);
		}
	}


	/**
	 * Returns a copy of the best individual in the population
	 */
	private Individual getBest() {
		best = null;
		for (Individual individual : population) {
			if (best == null) {
				best = individual.copy();
			} else if (individual.fitness < best.fitness) {
				best = individual.copy();
			}
		}
		return best;
	}

	/**
	 * INITIALISATION. Generates a randomly initialised population
	 */
	private ArrayList<Individual> initialise() {
		population = new ArrayList<>();
		for (int i = 0; i < Parameters.populationSize; ++i) {
			// chromosome weights are initialised randomly in the constructor
			Individual individual = new Individual();
			population.add(individual);
		}
		evaluateIndividuals(population);
		return population;
	}
	
	
	private ArrayList<Individual> PosNegInitialise() {
		population = new ArrayList<>();
		for (int i = 0; i < Parameters.populationSize; ++i) {
			// chromosome weights are initialised randomly in the constructor
			Individual individual = new Individual();
			Individual individual2 = individual.copy();
			
			for (int j=0; j<individual2.chromosome.length; j++) {
				// Flip chromes
				individual2.chromosome[j] = 0 - individual2.chromosome[j];
			}
			individual.fitness = Fitness.evaluate(individual, this);
			individual2.fitness = Fitness.evaluate(individual2, this);
			
			if (individual.fitness < individual2.fitness) {
				population.add(individual);
			} else {
				population.add(individual2);
			}
		}
		return population;
	}
	
	
	
	/**
	 * SELECTION
	 */
	private Individual randSelect() {
		Individual parent = population.get(Parameters.random.nextInt(Parameters.populationSize));
		return parent.copy();
	}
	private Individual tournamentSelect() {
		/**
		 * Elitism - copy the best chromosome (or a few best chromosomes) to new population
		 * (happens if tournament size is equal to total pop size)
		 * 1 - Pick t solutions completely at random from the population
		 * 2 - Select the best of the t solutions to be a parent
		 */
		final int TOURNAMET_POPULATION_SIZE = Parameters.tournamentSize;
		
		Collections.shuffle(population);
		Individual parent = population
				.stream()
				.limit(TOURNAMET_POPULATION_SIZE)
				.sorted((c1, c2) -> c1.compareTo(c2))
				.findFirst()
				.orElse(null);
		return parent;
	}
	// Fitness proportionate selection - roulette wheel selection
	private Individual rouletteSelect() {
		// calculate the total weight
		double weightSum = population
				.stream()
				.mapToDouble(c -> 1 - c.fitness)
				.sum();
		
		// Generate a random number between 0 and weightSum
		double random_seed = weightSum * Parameters.random.nextDouble();
		// Find random value based on weights
		for(Individual indiv : population) {		
			random_seed -= (1 - indiv.fitness);		
			if(random_seed < 0) 
				return indiv;
		}
		return population.get(-1);
	}

	
	private Individual rankSelect() {
		double[] fitness = new double[Parameters.populationSize];
		for (int i = 0; i < Parameters.populationSize; i++) {
	        fitness[i] = i + 1;
	    }

	    unitize1(fitness);
	    return population.get(random(fitness));
	}
	

    public static double norm1(double[] x) {
        double normalization = 0.0;
        for (double normalization_value : x) {
        	normalization += Math.abs(normalization_value);
        }
        return normalization;
    }
    
    public static void unitize1(double[] x) {
        double n = norm1(x);
        for (int i = 0; i < x.length; i++) {
            x[i] /= n;
        }
    }
    
    public static int random(double[] prob) {
        int[] ans = random(prob, 1);
        return ans[0];
    }

    public static int[] random(double[] prob, int n) {
        // set up alias table
        double[] q = new double[prob.length];
        for (int i = 0; i < prob.length; i++) {
            q[i] = prob[i] * prob.length;
        }

        // initialize a with indices
        int[] a = new int[prob.length];
        for (int i = 0; i < prob.length; i++) {
            a[i] = i;
        }

        // set up H and L
        int[] HL = new int[prob.length];
        int head = 0;
        int tail = prob.length - 1;
        for (int i = 0; i < prob.length; i++) {
            if (q[i] >= 1.0) {
                HL[head++] = i;
            } else {
                HL[tail--] = i;
            }
        }

        while (head != 0 && tail != prob.length - 1) {
            int j = HL[tail + 1];
            int k = HL[head - 1];
            a[j] = k;
            q[k] += q[j] - 1;
            tail++;                                  // remove j from L
            if (q[k] < 1.0) {
                HL[tail--] = k;                      // add k to L
                head--;                              // remove k
            }
        }

        // generate sample
        int[] samples = new int[n];
        for (int i = 0; i < n; i++) {
            double rU = Parameters.random.nextDouble() * prob.length;

            int k = (int) (rU);
            rU -= k;  

            if (rU < q[k]) {
            	samples[i] = k;
            } else {
            	samples[i] = a[k];
            }
        } 

        return samples;
        
    }

	private ArrayList<Individual> arithmeticCrossover(Individual first_parent, Individual second_parent){
		Individual child = new Individual();
		for (int i = 0; i < first_parent.chromosome.length; i++){
			double avgChrom = (first_parent.chromosome[i] + second_parent.chromosome[i]) / 2;
			child.chromosome[i] = avgChrom;
		}
		ArrayList<Individual> children_set = new ArrayList<>();
		children_set.add(child);
		return children_set;
	}
	
	/**
	 * REPLACEMENT
	 */
	private void replaceWorst(ArrayList<Individual> individuals) {
		for(Individual individual : individuals) {
			int idx = getWorstIndex();		
			population.set(idx, individual);
		}		
	}
	// Replace using tournament - same as selection but with worst
	private void tournamentReplace(ArrayList<Individual> individuals) {
		final int TOURNAMET_SIZE = Parameters.tournamentSize;
		
		for (Individual individual : individuals) {
			Collections.shuffle(population);
			Individual worstChrom = population
				.stream()
				.limit(TOURNAMET_SIZE)
				.sorted((c1, c2) -> c2.compareTo(c1))
				.findFirst()
				.orElse(null);

			population.remove(worstChrom);
			population.add(individual);
		}
	}
	
	
	private void constrainedMutation(ArrayList<Individual> individuals) {		
		for(Individual individual : individuals) {
			for (int i = 0; i < individual.chromosome.length; i++) {
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						double oldFitness = individual.fitness;
						individual.chromosome[i] += (Parameters.mutateChange);
						individual.fitness = Fitness.evaluate(individual, this);
						if (individual.fitness > oldFitness) {
							// revert if bad choice was made
							individual.chromosome[i] -= (Parameters.mutateChange);
						}
					} else {
						double oldFitness = individual.fitness;
						individual.chromosome[i] -= (Parameters.mutateChange);
						individual.fitness = Fitness.evaluate(individual, this);
						if (individual.fitness > oldFitness) {
							// revert if bad choice was made
							individual.chromosome[i] += (Parameters.mutateChange);
						}
					}
				}
			}
		}		
	}

	/**
	 * CROSSOVER 
	 */
	private ArrayList<Individual> uniformCrossover(Individual first_parent, Individual second_parent){
		Individual first_children = new Individual();
		Individual second_children = new Individual();
		
		for (int i = 0; i < first_parent.chromosome.length; i++){
			if(Parameters.random.nextBoolean()){
			first_children.chromosome[i] = first_parent.chromosome[i];
			second_children.chromosome[i] = second_parent.chromosome[i];
			} else {
		    first_children.chromosome[i] = second_parent.chromosome[i];
		    second_children.chromosome[i] = first_parent.chromosome[i];
			}
		}
		
		ArrayList<Individual> children_set = new ArrayList<>();
		children_set.add(first_children);
		children_set.add(second_children);	
		return children_set;
	}
	private ArrayList<Individual> onePointCrossover(Individual first_parent, Individual second_parent){
		Individual first_children = new Individual();
		Individual second_children = new Individual();
		int cutPoint = Parameters.random.nextInt(first_parent.chromosome.length);
		
		for (int i = 0; i < first_parent.chromosome.length; i++){
			if(i < cutPoint){
			first_children.chromosome[i] = first_parent.chromosome[i];
			second_children.chromosome[i] = second_parent.chromosome[i];
			} else {
		    first_children.chromosome[i] = second_parent.chromosome[i];
		    second_children.chromosome[i] = first_parent.chromosome[i];
			}
		}
		
		ArrayList<Individual> children_set = new ArrayList<>();
		children_set.add(first_children);
		children_set.add(second_children);	
		return children_set;
	}
	
	
	
	private ArrayList<Individual> twoPointCrossover(Individual first_parent, Individual second_parent){
		Individual first_children = new Individual();
		Individual second_children = new Individual();
		
		int chromosomeLength = first_parent.chromosome.length;
		int cutPoint1 = Parameters.random.nextInt(chromosomeLength);
		int cutPoint2 = Parameters.random.nextInt((chromosomeLength - cutPoint1) + 1) + cutPoint1;
		
		for (int i = 0; i < chromosomeLength; i++){
			if(i < cutPoint1 || i >= cutPoint2){
			first_children.chromosome[i] = first_parent.chromosome[i];
			second_children.chromosome[i] = second_parent.chromosome[i];
			} else {
		    first_children.chromosome[i] = second_parent.chromosome[i];
		    second_children.chromosome[i] = first_parent.chromosome[i];
			}
		}
		
		ArrayList<Individual> children_set = new ArrayList<>();
		children_set.add(first_children);
		children_set.add(second_children);	
		return children_set;
	}
	
	
	/**
	 *  Standard MUTATION
	 */
	private void mutate(ArrayList<Individual> individuals) {		
		for(Individual individual : individuals) {
			for (int i = 0; i < individual.chromosome.length; i++) {
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[i] += (Parameters.mutateChange);
					} else {
						individual.chromosome[i] -= (Parameters.mutateChange);
					}
				}
			}
		}		
	}

	/*
	 * getWorstIndex(): iterate the population list, and returns (if this exists, 
	 * 					the index of the worst performing solutions (population member); 
	 * @return: Index of the worst member found;  
	 */
	private int getWorstIndex() {
		Individual worst_member = null, current_population_member = null;
		int index_worst_memer_found = -1;
		for (int member_location = 0; member_location < population.size(); member_location++) {
			current_population_member = population.get(member_location);
			if (worst_member == null) {
				worst_member = current_population_member;
				index_worst_memer_found = member_location;
			} else if (current_population_member.fitness > worst_member.fitness) {
				worst_member = current_population_member;
				index_worst_memer_found = member_location; 
			}
		}
		return index_worst_memer_found;
	}
	

	@Override
	public double activationFunction(double x) {
		switch(Parameters.activationType) {
		case RELU:
			return (x > 0) ?  x : -1;
		case TANH:
			if (x < -20.0) {
				return -1.0;
			} else if (x > 20.0) {
				return 1.0;
			}
			return Math.tanh(x);
		default:
			return (x > 0) ? x : (Math.pow(Math.E, x) - 1) / 10;
		}
	}
}
