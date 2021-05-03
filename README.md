# Path Planning using Parallel Computing

The demand for Door-to-Door (D2D) delivery services has seen a huge surge in the wake of the pandemic and lockdowns. Path planning algorithms will play a huge role in the future when delivery is performed autonomously using intelligent UGVs and electric cars. We investigate a fast path planning technique that focuses on the core of Dijkstra’s algorithm and we implement the same on C with the augmentation of parallel computing platforms such as OpenMP and/or MPI. We also look into how such a technique can scale up for multiple agents. We simulate a small-scale path-planning scenario to visually demonstrate the algorithm. We also look into two different practical possibilities:
- The edges of the graph may be constrained to contain a particular number of agents beyond which collisions are guaranteed. One course of action is where we develop a multi-agent path planning approach to achieve trajectories that optimize the total time taken for all agents in such a scenario.
- While some agents are capable of moving to any given node due to their small physical requirements, others may not be able to find a path to the final location owing to larger restrictions. This provides the second course of action where we develop a multi-agent path planning approach to achieve optimal trajectories in such scenarios.

## Authors

* **Aaron John Sabu** - [aaronjohnsabu1999](https://github.com/aaronjohnsabu1999)
* **Athul CD** - [athulcd](https://github.com/athulcd)
* **Ayan Sharma** - [ayan-2004](https://github.com/ayan-2004)
* **Vaibhav Malviya** - [Vaibhav110](https://github.com/Vaibhav110)

<p align='center'>Created with :heart: by ...us...</p>
