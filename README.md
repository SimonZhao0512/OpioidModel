# OpioidModel
Program for modeling the prescription opioid epidemic for given data:  **APAM SEM 4903 Columbia University**

Opioid Model Documentation:
by Simon,Roland and Noah
10.16.19

### 1. **HOW TO USE:**
a. Make sure OpioidModel file contains the following files:
	opioidmodelmain.py
	opioidmodelmainnew.R
	opioidmodule.py
	properly formatted csv file  
b. Make sure that both opioidmodel.txt and the OpioidModel folder 
   are in your main   directory (i.e. ~/user/) OR change the output/input directories directly in both the scripts
for R and python  
c. From your Terminal, run ./opioidmodel.txt
	

**NOTE:** For a properly formatted csv file, please refer to the
      inputformat csv file in OpioidModel.
	

### 2. **OUTPUT:**
Based on your selections, a OpioidModelResult.pdf will be created in the  OpioidModel folder, containing your desired plot.

### 3. **OPERATION:**  
The program can perform the following computations and output a plot or/and an csv dataset:  
  a. Prescription opioid overdose attributed to the addicted class per year    
  b. Parameter(1) correlation plot for combinations resulting in OD deaths that were statistically indifferent(2)
     from the true od deaths in 2017  
  c. Prediction of 2017 od deaths using best parameters(3) for 2011-2016  
  d. Sobol Numerical Sensitivity Analysis(4): Finding first order and total order sensitivities of populations
     to our parameters  
  
  **NOTES:**  
  1. Parameters varied for options a-c: Prescription rate (S -> P), 
   rate of ending prescription without addiction (S -> P) and treatment entry rate (A -> R)  
  2. Results failed to reject a null hypothesis of no difference in a student t-test for a minimum sample of 50 values
     and a significance level of p < 0.1  
  3. Best parameters defined as those resulting in the minimal difference from the true od deaths for given year  
  4. This procedure uses sobol and saltelli modules from [SALib](https://salib.readthedocs.io/en/latest/api.html#sobol-sensitivity-analysis)
	

For any other questions, please contact Noah @ nmi2106@columbia.edu


