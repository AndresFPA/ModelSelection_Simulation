Code used for the simulation study in the paper "Selecting the number of clusters in Mixture Multigroup Structural Equation Modeling". A quick description of each file is listed below:

# Simulation Functions
**DataGeneration**: Script for the data generation.

**do_sim**: main script to run the simulation.

## Analyses & Evaluation Files
The analysis files process the results to bring a summary of the data in the form of a table or figure.

**Analysis**: Main analysis of the simulation. Contains the main results regarding the performance of the six model selection measures.

**Analysis2ndBest**: Analysis to determine the model selection performance including the best two models for all measures.

**Analysis_2to6**: Re-analysis of all model selection measures (excluding CHull) when excluding the 1-cluster model.

**AnalysisFacLogL**: Model selection performance using an alternative loglikelihood (results in the paper's supplemental material).

**AnalysisIgn**: Analysis when ignoring non-invariances. Not included in the paper.

The evaluation files evaluate the performance of the models (e.g., did we under- or over-selected?).

The corresponding evaluation files have the same names (e.g., evaluation2ndBest).

## Functions
It contains the functions required for the simulation, such as the main function of MMG-SEM or the ModelSelection function.

## Entropy
Contains the prior entropy analysis used to know the cluster separability of the generated data in the simulation.

## Results 
Contains the results of each one of the conditions evaluated in the simulation.
