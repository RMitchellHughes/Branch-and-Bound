# Branch-and-Bound-(Tensorflow-Implementation)
This repository provides open source software to support the methodologies of the paper "Branch and Bound to Assess Stability of Regression Coefficients in Uncertain Models" by Knaeble, Hughes, Rudolph, Abramson, and Razo (see https://arxiv.org/abs/2408.09634).

The current branch (Branch-and-Bound-(Tensorflow-Implementation)) is intented to showcase a Tensorflow adaptation to my original code (found in the main branch). More than anything, the Tensorflow version of the algorithm is intended to be a proof of concept rather than a heavily optimized algorithm. The main changes were the implementation of a linear regression model (utilizing conjugate gradient methods for least squares parameter estimation) and alterations to my original branch and bound code to to make use of that model.

**Disclaimer:** This is intended to be a seperate yet complimentary project to the original paper found above. More information about the original algorithm and dataset used can be found in the paper.

**File Details:**

The dataset is given in the **"NHANES_07_12.csv"** file. The example application uses this dataset.

In the file **"BB_confound_tf.py"**, two classes are found:
- _Node class_, represents an individual model that our algorithm checks.
- _BB_confound class_, takes explanatory variable x, response variable y, and a matrix of covariates s as inputs. x and y should be pandas dataframes with one column. s is a pandas dataframe also with multiple columns. Returns the maximum and minimum slope coefficients for x over the space of all possible models. Use the print function for more details about runtimes, nodes checked, which models resulted in the max/min, etc.

In the file **"LinearRegressionModel.py"**, one class is found:
- _LinearRegressionModel class_, represents a linear regression model. Takes an mxn tensorflow tensor of features (the intercept column of ones is added automatically) and a mx1 tensorflow tensor of the target variable. Will estimate least squares parameters via the conjugate gradient method, calculate residuals, and calculate coefficient of determination upon request (parameter estimation is the only thing done automatically).

One more class is required to run the code:
- A final class called _ConfoundingInterval_ written by Mark Abramson and supplied by Brian Knaeble is used to compute the confounding intervals. To use the method from this class, the file ConfoundingInterval.py should be in the working directory. It can be downloaded from https://github.com/bknaeble/ConfoundingIntervals/tree/master.

An example of how the function works is given in the file **example_tf.py**. It assumes that the files BB_confound_tf.py, LinearRegressionModel.py, ConfoundingInterval.py and NHANES_07_12.csv are found in the working directory.
