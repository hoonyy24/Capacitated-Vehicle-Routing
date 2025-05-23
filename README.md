# Capacitated-Vehicle-Routing

In this project, I worked with XML100 data, where the optimal solutions were computed using a Vehicle Routing Problem (VRP) solver. I developed a Capacitated Vehicle Routing Problem (CVRP) solution based on the Adaptive Large Neighborhood Search (ALNS) approach. To get started, I used the Clarke–Wright Savings algorithm to create an initial solution. This involved assigning each customer a direct route from the depot and then merging routes based on savings values, all while keeping vehicle capacity constraints in check.

Next, I improved the solution through an iterative process. In the ALNS loop, I applied a combination of removal and insertion operators, along with local search techniques, to enhance the results. I selected operators using a Softmax-based method and incorporated a probabilistic acceptance strategy—guided by a temperature parameter T—to occasionally accept slightly worse solutions, helping me escape local optima. When improvements stalled, I introduced a diversification strategy to jump to new starting points. Finally, I fine-tuned the solution using Variable Neighborhood Descent (VND) and Variable Neighborhood Search (VNS) techniques to get as close as possible to the optimal outcome.

For evaluation, I tested my approach on 1,000 randomly selected problems from a total of 10,000 in the XML100 instance from CVRPLIB. The VRP solver provided the benchmark optimal solutions, and my calculations came in with an average cost deviation of just 4.65% from those optimal results.


## An example of an optimal route computed by the ALNS algorithm
![RoutePlot_XML100_2376_15_vrp](https://github.com/user-attachments/assets/fbdcabc3-90b4-4df2-911e-e03dae563918)
