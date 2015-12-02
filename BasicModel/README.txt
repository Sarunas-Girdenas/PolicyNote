In this folder there are files that were used to run basic model (with different Taylor type rules). 

The files are:

welfareLoopEta - function that takes eta (bargaining power) parameter along with simulation results from dynare and 
		 a struct of parameters values. Returns welfare.
parameterDefinition - script that is used to define parameter values
steadyStateFunction - function that takes parameters and returns steady state values of variables
extendModelUnLoop_steadystate - function that supplies steady state solution to dynare file
mainSimulFileUn - script full of tiny and useful code snippets that automate some simulation tasks.
welfareLoop - the same function as welfareLoopEta but does not take eta as an input
extendModelUnLoop - dynare file of the model

Sarunas Girdenas, University of Exeter, UK, sg325@exeter.ac.uk

