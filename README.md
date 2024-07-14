# Branch-and-Bound
This repository provides open source software to support the methodologies of the paper "Branch and Bound to Assess Stability of Regression Coefficients in Uncertain Models" by Knaeble, Razo, Hughes, and Abramson

In the file "BB_code.R", four functions are found:
- calc.I.res, this function calculates the matrix of residuals from regressing each vector in I onto J.
- f, this function computes the confounding intervals. It was supplied by Knaeble with more details at https://github.com/bknaeble/ConfoundingIntervals/tree/master.
- BB.confound, this function is our branch and bound (BB) algorithm. It takes a vector x, a vector y, and a matrix w as inputs. It outputs the maximum and minimum slope coefficients for x over the space of all possible models.
- BB.confound.reorder, this function is our branch and bound with reordering (BBR) algorithm. It also takes a vector x, a vector y, and a matrix w as inputs. It outputs the maximum and minimum slope coefficients for x over the space of all possible models.

In the file "NHANES_07_12.csv", the data used for section 4 of the paper is provided. Information about how the data was wrangled and cleaned is given in both section 4 and Appendix B.

An example of how the function works is given in the file "example_application.R". It assumes that the functions from BB_code.R have been loaded and the file NHANES_07_12.csv is found in the working directory.
