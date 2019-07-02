Codes used for the main analysis are stored here.

1. **Cleaning.R**: Here, the data is prepared for the analysis. Beside cleanings specific to the data, exposure variables are firstly transformed using Box-Cox transformation. Next, Z-transformation has been done to make the effect sizes comparable.

2. **EWAS.R**: the data cleaned in the previous step is used here to run an EWAS. Functiom *cont_ewas* is the main function, where Cox model is used to associate an exposure with time to failure (CHD), controlling for adjusitng varaibles. This function checks for proportional hazard assumption as well. Also, it calculates the statistical power of the association tested. Finally, this function calculates the VIF to check the presence of multi-collinearity. The function outputs the results from Cox model, p-value related to proportional hazard assumption, statistical power, and the VIF.

3. **Permutation.R**: 
