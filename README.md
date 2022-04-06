# Rocket landing
The current repository contains the SET10107 Computational intelligence coursework, given by Edinburgh Napier University.

## Technologies
![Java](https://img.shields.io/badge/java-%23ED8B00.svg?style=for-the-badge&logo=java&logoColor=white)

## Goals
* Prune training parameters to optimize the weights of a Multi-layer Multi-layer Perceptron, applied to the landing system of a rocket 
* Get lowest fitness score, produced by the average of multiple runs of the same Evolutionary algorithm; 

##  Notes:
Actions has been taken on the following files:
* ExampleEvolutionaryAlgorithm.java
* Parameters.java

## Evaluation
When running StartGui.java, the running of the UI will slow down the execution of the algorithm, reporting the trends of fitness score just until the 10.000th running;
Instead, this project runs the EA 20.000 times; To optimize the training time, and visualize final results, you may want to start the training by running StartNoGui.java
