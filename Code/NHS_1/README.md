Codes used for the main analysis are stored here.

1. **Cleaning.R**: Here, the data is prepared for the analysis. Beside cleanings specific to the data, exposure variables are firstly transformed using Box-Cox transformation. Next, Z-transformation has been done to make the effect sizes comparable.

2. **EWAS.R**: The data cleaned in the previous step is used here to run an EWAS. Functiom *cont_ewas* is the main function, where Cox model is used to associate an exposure with time to failure (CHD), controlling for adjusitng varaibles. This function checks for proportional hazard assumption as well. Also, it calculates the statistical power of the association tested. Finally, this function calculates the VIF to check the presence of multi-collinearity. The function outputs the results from Cox model, p-value related to proportional hazard assumption, statistical power, and the VIF.

3. **Permutation.R**: Here, using adjusting varibales, the probability of developing CHD at each time for each participant is calculated, using Cox model. The probabilities are used to shuffle the CHD status. After each shuffling, the association between each exposure and time to failure (CHD) is recalculated and the related p-value is collected. The shuffling is repeated 1,000 times. Finally, 1,0000 null p-values are collected for each exposure.

4. **FDR_estimation.py**: The outcome of the previous step is used here to estimate the FDR, to control for type I error due to multiple hypotheses testing. The FDR is estimated to be the ratio of the proportion of results that were called significant at a given level $\alpha$ in the null distribution and the proportion of results called significant from real tests.
