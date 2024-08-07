# Branch-and-Bound
This repository provides open source software to support the methodologies of the paper "Branch and Bound to Assess Stability of Regression Coefficients in Uncertain Models" by Knaeble, Hughes, Rudolph, Abramson, and Razo

In the file "BB_code.R", four functions are found:
- calc.z.res, this function returns the matrix of residuals from regressing each vector in z onto the span of w.
- f, this function computes the confounding intervals. It was supplied by Knaeble with more details at https://github.com/bknaeble/ConfoundingIntervals/tree/master.
- BB.confound, this function is our branch and bound algorithm. It also takes an explanatory vector x, a response vector y, and a matrix of covariates s as inputs. It outputs the maximum and minimum slope coefficients for x over the space of all possible models.

In the file "NHANES_07_12.csv", the data used for section 4 of the paper is provided. Information about how the data was wrangled and cleaned is given in both section 4 and Appendix B of the paper.

An example of how the function works is given in the file "example_application.R". It assumes that the functions from BB_code.R have been loaded and the file NHANES_07_12.csv is found in the working directory.
