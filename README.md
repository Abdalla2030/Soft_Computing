# Soft_Computing

Assignments of Soft Computing/Genetic Algorithms Course ( CS465 )

# Assignments
***************************************************************************************************
### Assignment-1 ( Knapsack Problem )

* About the problem:
The knapsack problem is a very well-known optimization problem. Given a knapsack that can carry weights up to a certain amount, and a number of items; each item has a weight and a value, we want to select the items to carry in the knapsack in order to maximize the total value.

* What you are required to do:
Write a genetic algorithm to solve the knapsack problem.

* What the input looks like: You’ll be given an input file with the following format:
• First line: Number of test cases (must be at least 1)
For each test case:
• Size of the knapsack
• Number of items
For each item:
• Weight and value separated by space
## Example 
```
2
10
3
4 4
6 6
5 3
15
5
12 4
1 2
4 10
1 1
2 2
```
* Important remarks to help you solve the problem:
1. Use a binary, one-dimensional chromosome.
2. The population size and initialization method you use are up to you. You can actually try different population sizes to see how this will affect your results. The
maximum number of generations is also up to you.
3. Think about how you will handle infeasible solutions. Infeasible solutions are solutions that violate the constraints of the problem; therefore, they are not
allowed.
4. Use rank selection and one-point crossover. Choose the methods of mutation and replacement that you find appropriate.
5. The output should consist of the test case index, the number of selected items, the total value, and the weight and value of each selected item
***************************************************************************************************
### Assignment-2 ( Smooth Curve Fitting Problem )

* About the problem:
Curve fitting is the process of constructing a curve, or mathematical function (polynomial
equation) that has the best fit to a series of data points, possibly subject to constraints. In smooth curve fitting, the function is constructed to
approximately fit the data. Given a set of points, we would like to fit a curve to them using a polynomial equation. 

* What you are required to do:
Write a genetic algorithm to find the best coefficients that would make the distance between the polynomial function and the points minimum.

* What the input looks like:
You’ll be given an input file with the following format:
• First line: Number of datasets (must be at least 1)
For each dataset:
• Number of data points and polynomial degree separated by space
For each point:
• x-coordinate and y-coordinate separated by space
## Example 
```
1
4 2
1 5
2 8
3 13
4 20
```
Note: In the example above, we have 1 dataset containing 4 points and we want to fit a 2nd
degree polynomial (a0, a1, a2).

* Important remarks:
1. The first thing you need to do is to think about how you will encode the solution in your chromosome and what the objective function will be.
2. You should use a floating-point representation of the chromosome.
3. You can try different population sizes to see how this will affect your results. The maximum number of generations is up to you.
4. Initialize the genes such that their values are in the range [-10,10].
5. Implement tournament selection, 2-point crossover, non-uniform mutation, and elitist replacement.
6. You must read the input from the given file and write the output to a file. The output should consist of the dataset index, the coefficients of the polynomial function, and its mean square error.
