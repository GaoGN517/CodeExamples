This is actually a special case of pairwise comparisons for post-ANOVA analyses, 
so basically, Bonferroni procedure works. To do pairwise comparisons, we have 
Scheffe's method, Fisher's LSD, Bonferroni method in general. 

Fisher's LSD: performs F test, if rejects the null, then perform t-test for all 
pairs each at alpha level. Not a good method for multi-comparisons.

Bonferroni: correction for the significance level - alpha / (number of comparisons). 

Bonferroni vs Fisher: Fisher controls the alpha-level error rate for each pairwise
comparison, but not the family error rate. Bonferroni does not control the family
error rate, but reduce the significance level based on the number of comparisons.
So Bonferroni confidence intervals are wider than Fisher's.

Scheffe's method: investigate all possible contrasts of the means. If F-test rejects
the null hypothesis at level alpha, then at least one contrast should be rejected in
this test. 

Tukey's Studentized method divide all pairs of means by the estimated standard
deviation of the mean, and compare them with tabled critical values. It is 
slightly better than Bonferroni in all-pairwise comparisons. Tukey's procedure 
is exact for equal sample sizes. The approximate procedure for unequal sizes is 
called Tukey-Kramer. 

Dunnett’s procedure is designed to compare each treatment to a control. For n 
groups, there are n - 1 comparisons. It is a bit more conservative than Bonferroni,
and it is an exact procedure.












