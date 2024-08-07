# Branch-and-Bound
This repository provides open source software to support the methodologies of the paper "Branch and Bound to Assess Stability of Regression Coefficients in Uncertain Models" by Knaeble, Hughes, Rudolph, Abramson, and Razo.

In the file **"NHANES_07_12.csv"**, the data used for section 4 of the paper is provided. Information about how the data was wrangled and cleaned is given in both section 4 and Appendix B of the paper.

**R Code:**

In the file **"BB_code.R"**, four functions are found:
- _calc.z.res_, this function returns the matrix of residuals from regressing each vector in z onto the span of w.
- _f_, this function computes the confounding intervals. It was written by Mark Abramson and supplied by Brian Knaeble with more details at https://github.com/bknaeble/ConfoundingIntervals/tree/master.
- _BB.confound_, this function is our branch and bound algorithm. It also takes an explanatory vector x, a response vector y, and a matrix of covariates s as inputs. It returns the maximum and minimum slope coefficients for x over the space of all possible models.

An example of how the function works is given in the file **"example.R"**. It assumes that the functions from BB_code.R have been loaded and the file NHANES_07_12.csv is found in the working directory.

**Python Code:**

In the file **"BB_confound.py"**, two classes are found:
- _Node class_, represents an individual model that our algorithm checks. It has parameters to store residuals of x, y, and z, as well as the subindices I_w and I_z.
- _BB_confound class_, takes explanatory variable x, response variable y, and a matrix of covariates s as inputs. x and y should be pandas dataframes with one column. s is a pandas dataframe also. Returns the maximum and minimum slope coefficients for x over the space of all possible models. Use the print function for more details about runtimes, nodes checked, which models resulted in the max/min, etc.
- A third class called _ConfoundingInterval_ written by Mark Abramson and supplied by Brian Knaeble is used to compute the confounding intervals. To use this function, the file ConfoundingInterval.py should be in the working directory. It can be downloaded from https://github.com/bknaeble/ConfoundingIntervals/tree/master.

An example of how the function works is given in the file **example.py**. It assumes that the files BB_confound.py, ConfoundingInterval.py and NHANES_07_12.csv are found in the working directory.
