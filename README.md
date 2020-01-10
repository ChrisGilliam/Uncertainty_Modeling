## Uncertainty_Modeling
Matlab code and functions for the testing scenarios analysed in "[A tutorial on uncertainty modeling for machine reasoning](https://doi.org/10.1016/j.inffus.2019.08.001)".

Increasingly we rely on machine intelligence for reasoning and decision making under uncertainty. The tutorial reviews the prevalent methods for model-based autonomous decision making based on observations and prior knowledge, primarily in the context of classification. Both observations and the knowledge-base available for reasoning are treated as being uncertain. Accordingly, the central themes of this tutorial are quantitative modeling of uncertainty, the rules required to combine such uncertain information, and the task of decision making under uncertainty. The paper covers the main approaches to uncertain knowledge representation and reasoning, in particular:
- Bayesian probability theory
- Possibility theory
- Dempster-Shafer belief functions theory (using transferrable belief model)
- Imprecise probability theory 

These approaches are illustrated on several testing scenarios as outlined below:

Testing Scenario 1 - Simple Classification
- _Script_1_ implements a Bayesian classifier
- _Script_5_ implements a possibilistic classifier
- _Script_9_ implements a transferrable belief model classifier

Testing Scenario 2 - Classification under Randomness via Monte-Carlo Simulations
- _Script_2_ implements a Bayesian classifier
- _Script_6_ implements a possibilistic classifier
- _Script_10_ implements a transferrable belief model classifier

Testing Scenario 3 - Imprecise Likelihood Specification
- _Script_3_ implements a standard Bayesian approach
- _Script_4_ implements Mahler's approach using random sets
- _Script_7_ implements a possibilistic solution
- _Script_11_ implements a transferrable belief model solution
- _Script_14_ implements a imprecise probabilistic solution

Testing Scenario 4 - Model Mismatch
- _Script_8_ compares the Bayesian and possibilistic classifiers when using an incorrect confusion matrix (Monte-Carlo runs)
- _Script_12_ compares the Bayesian, possibilistic and transferrable belief model classifiers when using an incorrect confusion matrix (Monte-Carlo runs)
- _Script_15_ solves a modified version of the scenerio using Imprecise probability theory - rather than an incorrect, the confusion matrix is imprecise.

Testing Scenario 5 - Dealing with different subspaces 
- _Script_13_ implements a transferrable belief model solution

## Authors:
Scripts 1 to 13 and 15 were written by B. Ristic.
Script 14 was written by A. Benavoli

## References
1) B. Ristic, C. Gilliam, M. Byrne, and A. Benavoli, "[A tutorial on uncertainty modeling for machine reasoning](https://doi.org/10.1016/j.inffus.2019.08.001)", Information Fusion, Vol 55, pp. 30-44, 2020.

