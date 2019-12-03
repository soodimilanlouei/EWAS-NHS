<p align="justify"> In *data_sample.sas7bdat*, a random sample is stored to show the format of the data that we used in our analysis. Please note that this is not a subset of the real data but an expression of how the real data may look like and what the required data format is for this analysis.
</p>


**Codes used for the main analysis are stored here.**

1. **Cleaning.R** <p align="justify"> Here, the data is prepared for the analysis. Beside cleanings and manipulation specific to the data, exposure variables are firstly transformed using Box-Cox transformation. Next, Z-transformation has been done to make the effect sizes comparable.
</p>

2. **EWAS.R** <p align="justify"> The data cleaned in the previous step is used here to run an EWAS. Functiom *cont_ewas* is the main function, where Cox model is used to associate an exposure with time to failure (CHD), controlling for adjusitng varaibles. This function checks for proportional hazard assumption as well. Also, it calculates the statistical power of the association tested. Finally, this function calculates the VIF to check the presence of multi-collinearity. The function outputs the results from Cox model, p-value related to proportional hazard assumption, statistical power, and the VIF. The expected output should be similar to Table 3 in the SI (without the FDR column).
</p>

3. **Permutation.R**  <p align="justify"> Here, using adjusting varibales, the probability of developing CHD at each time for each participant is calculated, using Cox model. The probabilities are used to shuffle the CHD status. After each shuffling, the association between each exposure and time to failure (CHD) is recalculated and the related p-value is collected. The shuffling is repeated 1,000 times. Finally, 1,000 null p-values are collected for each exposure. The output should look like below with 1000 rows.
</p>

<img src="https://github.com/soodimilanlouei/EWAS-NHS/blob/master/Code/NHS_1/permutation-output.png" width="700">

4. **FDR_estimation.py**  <p align="justify"> The outcome of the previous step is used here to estimate the False Discovery Rate (FDR), to control for type I error due to multiple hypotheses testing. The FDR is estimated to be the ratio of the proportion of results that were called significant at a given level alpha in the null distribution and the proportion of results called significant from real tests. The output should look like Table 3 in the SI with columns Exposure, P-value, and FDR.
</p>

5. **Correlation_Analysis.R**  <p align="justify"> The correlation analysis is done in the SI. The Spearman correlation between significant associations found in NHS. Similar to (Patel and Manrai, 2014), a permutation-based approach is used to estimate the two-sided P-value of significance for each pair of correlations. Given a pair of X and Y, we shuffled values of X and computed the correlation with Y, repeating this procedure for all pairs. The P-value for a correlation was the fraction of correlations from the permuted dataset with greater absolute value. The output should look like Figure 3 in the SI.
</p>

6. **Plots.py**  <p align="justify"> This file contains the code for creating most of the plots in the manuscript.
</p>

7. **main_SESfigures.m** <p align="justify"> The Socio-economics (SES) analysis is done in the SI. We analyzed the effect of SES and compared it with our primary results. This file contains the code for producing Figures 4 and 5 in the SI.
</p>
