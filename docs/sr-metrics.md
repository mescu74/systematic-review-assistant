# Calculating Accuracy Metrics for Systematic Review Screening: Formulas and Applications

Systematic reviews require rigorous screening of literature to identify relevant studies. Evaluating the performance of this screening process—whether conducted manually or with automated tools—relies on specific statistical metrics. This report details how these metrics are computed, their formulas, and their significance in systematic review methodology.

## The Foundation: Confusion Matrix

At the heart of screening evaluation is the confusion matrix, which categorizes screening decisions into four fundamental outcomes [1] [2] [3]:

|                          | Articles Actually Relevant | Articles Actually Irrelevant |
| :----------------------- | :------------------------- | :--------------------------- |
| **Predicted Relevant**   | True Positives (TP)        | False Positives (FP)         |
| **Predicted Irrelevant** | False Negatives (FN)       | True Negatives (TN)          |

This $2 \times 2$ classification forms the basis for calculating all performance metrics in systematic review screening. Each cell represents an important outcome that affects the review's reliability and efficiency.

## Primary Screening Performance Metrics

### Sensitivity and Specificity

Sensitivity (also called recall) measures the ability to identify truly relevant articles, which is particularly crucial in systematic reviews [4] [5]:

$$ \text{Sensitivity} = \frac{\text{TP}}{\text{TP} + \text{FN}} $$

Sensitivity represents the probability that a relevant article will be correctly identified during screening. In medical screening contexts, this is equivalent to the probability of a positive test result in a case with disease [14].

Specificity measures the ability to correctly identify irrelevant articles [4] [6]:

$$ \text{Specificity} = \frac{\text{TN}}{\text{TN} + \text{FP}} $$

Specificity indicates the probability that an irrelevant article will be correctly excluded during screening.

## Predictive Values

Positive Predictive Value (PPV), also known as precision, measures the proportion of articles identified as relevant that are truly relevant [7] [8]:

$$ \text{PPV} = \frac{\text{TP}}{\text{TP} + \text{FP}} $$

PPV is particularly important for efficiency, as it indicates how much of the screener's effort in reviewing "relevant" articles is well-directed [9].

Negative Predictive Value (NPV) measures the proportion of articles identified as irrelevant that are truly irrelevant [7] [10]:

$$ \text{NPV} = \frac{\text{TN}}{\text{TN} + \text{FN}} $$

A high NPV indicates minimal loss of relevant articles during screening.

### Accuracy

Accuracy represents the overall proportion of correct classifications [11]:

$$ \text{Accuracy} = \frac{\text{TP} + \text{TN}}{\text{TP} + \text{TN} + \text{FP} + \text{FN}} $$

While intuitive, accuracy can be misleading in systematic reviews due to class imbalance (typically many more irrelevant than relevant articles). Metrics like the F1 Score or Matthews Correlation Coefficient are often more informative in such scenarios.

## Advanced Performance Metrics

### F1 Score

The F1 score represents the harmonic mean of precision and recall, providing a balance between these two metrics, and is often more useful than accuracy for imbalanced datasets [1] [3].

$$ F1 = 2 \times \frac{\text{Precision} \times \text{Recall}}{\text{Precision} + \text{Recall}} $$

This metric is valuable when seeking a balance between finding all relevant articles (recall) and minimizing the screening of irrelevant ones (precision).

### Matthews Correlation Coefficient

The Matthews Correlation Coefficient (MCC) provides a balanced measure of screening quality, particularly effective even with highly imbalanced classes [11] [12]:

$$ \text{MCC} = \frac{\text{TP} \times \text{TN} - \text{FP} \times \text{FN}}{\sqrt{(\text{TP}+\text{FP})(\text{TP}+\text{FN})(\text{TN}+\text{FP})(\text{TN}+\text{FN})}} $$

MCC values range from -1 to +1, where +1 represents perfect prediction, 0 indicates random prediction, and -1 signifies inverse prediction. The MCC is particularly valuable for systematic reviews because it handles class imbalance well [3].

## Agreement Statistics

When evaluating the consistency between different screeners or between manual and automated screening, agreement statistics become important.

### Cohen's Kappa

Cohen's Kappa measures inter-rater reliability adjusted for chance agreement [13] [14]:

$$ \kappa = \frac{p_o - p_e}{1 - p_e} $$

Where:

- $p_o$ is the observed agreement proportion (i.e., $\frac{\text{TP} + \text{TN}}{\text{Total}}$)
- $p_e$ is the expected agreement proportion due to chance

### Prevalence and Bias Adjusted Kappa (PABAK)

Prevalence and Bias Adjusted Kappa (PABAK) adjusts for bias and prevalence effects, reporting what Kappa would be if prevalence were 50% and there were no bias [15] [3]:

$$ \text{PABAK} = 2p_o - 1 $$

Where $p_o$ is the observed agreement proportion (i.e., $\frac{\text{TP} + \text{TN}}{\text{Total}}$). This simplified formula is derived assuming a $2 \times 2$ table where expected agreement ($p_e$) is set to 0.5.

## Likelihood Ratios

Positive Likelihood Ratio (LR+) and Negative Likelihood Ratio (LR-) are particularly useful in diagnostic contexts for assessing the value of a screening test [10]:

$$ \text{LR+} = \frac{\text{Sensitivity}}{1 - \text{Specificity}} $$

$$ \text{LR-} = \frac{1 - \text{Sensitivity}}{\text{Specificity}} $$

## Application in Systematic Review Practice

### Prioritizing Sensitivity

In systematic review screening, sensitivity is generally prioritized over specificity or precision [15]. Research indicates that a recall (sensitivity) over 0.99 should be targeted to ensure comprehensive literature coverage [5]. This prioritization stems from the fundamental principle that missing relevant studies introduces more serious bias than including some irrelevant ones initially.

### Performance in Modern Automated Screening

With the emergence of machine learning and large language models for systematic review automation, these metrics have taken on additional importance for evaluating their performance [16] [3]. Advanced models like GPT-4 variants and Mistral models have demonstrated promising performance across multiple metrics, with Matthews Correlation Coefficients ranging from 0.290 to 0.342 and PABAK scores between 0.537 and 0.621 in some studies [3]. The application of Random Forest classifiers has shown that more flexible approaches to inclusion criteria can sometimes yield better screening performance than requiring all criteria to be met simultaneously [3].

## Conclusion

Accurate evaluation of systematic review screening performance requires appropriate application of these statistical metrics. While sensitivity remains paramount to ensure comprehensive coverage, other metrics like MCC, F1 Score, and PABAK offer robust evaluation, particularly in the context of class imbalance inherent in screening tasks. As automated screening tools continue to evolve, these metrics will play an increasingly important role in validating their performance and ensuring the reliability of systematic reviews.

For researchers conducting or evaluating systematic reviews, understanding these metrics and their interrelationships is essential for assessing screening quality and making informed methodological choices that balance comprehensiveness with resource efficiency.

## References

[1] [https://www.sciencedirect.com/topics/engineering/confusion-matrix](https://www.sciencedirect.com/topics/engineering/confusion-matrix)
[2] [https://www.v7labs.com/blog/confusion-matrix-guide](https://www.v7labs.com/blog/confusion-matrix-guide)
[3] [https://www.pnas.org/doi/10.1073/pnas.2411962122](https://www.pnas.org/doi/10.1073/pnas.2411962122)
[4] [https://methods.cochrane.org/sites/methods.cochrane.org.sdt/files/uploads/Chapter%2010%20-%20Version%201.0.pdf](https://methods.cochrane.org/sites/methods.cochrane.org.sdt/files/uploads/Chapter%2010%20-%20Version%201.0.pdf)
[5] [https://pmc.ncbi.nlm.nih.gov/articles/PMC9277646/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9277646/)
[6] [https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2017.00307/full](https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2017.00307/full)
[7] [https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values](https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values)
[8] [https://geekymedics.com/sensitivity-specificity-ppv-and-npv/](https://geekymedics.com/sensitivity-specificity-ppv-and-npv/)
[9] [https://pmc.ncbi.nlm.nih.gov/articles/PMC6664615/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6664615/)
[10] [https://www.ncbi.nlm.nih.gov/books/NBK557491/](https://www.ncbi.nlm.nih.gov/books/NBK557491/)
[11] [https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html)
[12] [https://www.omnicalculator.com/statistics/matthews-correlation-coefficient](https://www.omnicalculator.com/statistics/matthews-correlation-coefficient)
[13] [https://real-statistics.com/reliability/interrater-reliability/cohens-kappa/](https://real-statistics.com/reliability/interrater-reliability/cohens-kappa/)
[14] [https://en.wikipedia.org/wiki/Cohen%27s_kappa](https://en.wikipedia.org/wiki/Cohen%27s_kappa)
[15] [https://stats.stackexchange.com/questions/19601/adjusting-kappa-inter-rater-agreement-for-prevalence](https://stats.stackexchange.com/questions/19601/adjusting-kappa-inter-rater-agreement-for-prevalence)
[16] [https://www.medrxiv.org/content/10.1101/2024.06.03.24308405v1.full-text](https://www.medrxiv.org/content/10.1101/2024.06.03.24308405v1.full-text)
