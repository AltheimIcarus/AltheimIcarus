# l2Match Filter-and-Verification Subgraph Isomorphism Algorithm
## Introduction
We proposed l2Match algorithm for Non-Induced SI problem to optimize exisiting methods, specifically:
(1) Avoid scanning neighbors of a vertex with mismatching labels in the filtering step.
(2) Reduce redundant traversal in the Forward-Candidates-Generation-and-Backward-Candidates-Pruning (FCGBCP) filtering method.
(3) Reduce the expensive removal of invalid candidates in the CECI indexing data structure of the CECI algorithm during the filtering step.
(4) Reduce the exploration of redundant search branches during the enumeration step.

We evaluate the performance of l2Match against GraphQL, CFL, CECI, and DP-iso algorithms in a common
framework (Shixuan Sun and Qiong Luo. In-Memory Subgraph Matching: an In-depth Study. SIGMOD 2020) [[1]](https://github.com/RapidsAtHKUST/SubgraphMatching).

We report two types of analysis:

(1) average performance for all query that is solvable by every algorithm within the time limit;

(2) the number of solved queries for each algorithm when one or more algorithms exceeded time limit to solve any query.

Our paper is accepted by ICEIC2024. Kindly contact us should you have any enquiry. [Full Text Manuscript](https://github.com/AltheimIcarus/l2Match/blob/master/ChiQinCheng-ICEIC.pdf)

## Metrics
The metrics are as such:
(1) Time elapsed are measured in nanoseconds but reported in seconds unless specifically stated.

(2) Total candidate count is calculated as $\Sigma_{u\in V(Q)}|C(u)|$.

(3) The time limit is set to $0.3\times 10^{3}$ seconds or 5 minutes.

(4) The enumeration is stopped at $10^{6}$ embeddings (isomorphic subgraphs).

(5) As the number of edges are usually larger than the number of vertices, we evaluate the results over increasing number of query edges.

## Hardness vs Satisfiability
A query is satisfiable if there exist at least one subgraph in the data graph that is isomorphic to it, otherwise unsatisfiable.
The hardness of a query is defined as:
| Hardness | Description |
| :----------------: | :---------------------------------------: |
|l2MatchJR| l2Match algorithm optimized with novel Jump-and-Redo (JR) method |
| Easy (E) | if it can be solved by all algorithm within the time limit |
| Easy/Hard (EH) | if it can be solved by more than one but not all algorithm within the time limit |
| Hard (H) | if it can be solved by exactly one algorithm within the time limit |
| Unsolved (U) | if none of the algorithm can solve it within the time limit |

E, EH, and H queries are guaranted to be satisfiable. Yet, any U query may be satisfiable or unsatisfiable.


## Compile
Execute commands below to compile the source code in the source directory. This program can only be compiled in Linux environment as it contains SIMD codes that are native to specific OS.  

```zsh
mkdir build
cd build
cmake ..
make
```

## Test
You may edit the test.cpp file accordingly to design your own test script.


## Execute
Use the 'test.out' program in the 'build/matching' directory to execute any test cases.
You will have to supply 3 additional parameters:
(1) -d <path-to-data-graph>
(2) -q <path-to-query-graph>
(3) -filter <algorithm-to-run>
The maximum embeddings count and time limit can be set in (test.cpp):

```zsh
TimeOutException* timeout_e = new TimeOutException();
std::thread([&timeout_e]{
    std::this_thread::sleep_for(300s); // 5 minutes
    timeout_e->is_throwable = true;
}).detach();
```

and

```zsh
size_t output_limit = std::numeric_limits<size_t>::max();
output_limit = 1000000;
```

For example:
```zsh
./test.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter GQL
```

## Input
Both the input query graph and data graph are vertex-labeled.
Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted
as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. Note that we require that the vertex
id is started from 0 and the range is [0,N - 1] where V is the vertex set. The following
is an input sample. You can also find sample data sets and query sets under the test folder.

Example:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Algorithms included

Follwing algorithms are included in the test environment:

|Parameter of Command Line (-filter) | Description |
| :-----------------------------------: | :-------------: |
|l2MatchJR| l2Match algorithm optimized with novel Jump-and-Redo (JR) method |
|GQL| GraphQL algorithm |
|CFL| CFL-Match algorithm|
|DPiso| DP-iso algorithm optimized with Failing Set Pruning method |
|CECI| CECI algorithm |

More settings can be found in the defined macros in 'configuration/config.h'.


## Recommendation
Our method, l2Match achieves performance improvement over CECI by 19.22%, CFL by 13.11%, GQL by 45.21%, and DPiso by 37.79%.
Performance improvement is calculated as the percentage of decrease $p=100\cdot (t^{algo}-t^{l2Match})\div t^{algo}$ where $t$ is the average query time and $algo$ represents each competing algorithm.



## Datasets
The datasets used in our paper are gathered from two sources:

(1) In-Memory Subgraph Matching: an In-depth Study (Sun, S., Luo, Q., 2020): [datasets](https://hkustconnect-my.sharepoint.com/:u:/g/personal/ssunah_connect_ust_hk/EQnXTic0PK9Fo1gkdDZRKOIBFIyMeBTP5rbju2ZfQdj-QA?e=SfGa8X), specifically DBLP, Human, Patents, Yeast, and Youtube datasets.

(2) Experimental Evaluation of Subgraph Isomorphism Solvers (Solnon, C., 2019): [datasets](https://perso.liris.cnrs.fr/christine.solnon/SIP.html), specifically Erdős–Rényi graphs (ERP).

The file size of the datasets exceeds the upload limit of the GitHub repository. Hence, we do not include the datasets in this repository.

## Overall Performance
The average and standard deviation of query time over increasing query size $|V(Q)|$ are reported in the figures below. Notice that $|V(Q)|$ is used as the x-axis instead of the number of query edges $|E(Q)|$ as the spikes and troughs are densely compacted, which hinders viewability. The proposed algorithm, l2Match outperforms CECI by 19.22\%, CFL by 13.11\%, GQL by 45.21\%, and DPiso by 37.79\%. Performance improvement is calculated as the percentage of decrease $p=100\cdot (t^{algo}-t^{l2Match})\div t^{algo}$ where $t$ is the average query time and $algo$ represents each competing algorithm. In summary, the overall performance ranking is l2Match $>$ CFL $>$ CECI $>$ DPiso $>$ GQL.


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/Performance.png" width="500" alt="Avg. query time over |V(Q)|" />
    <figcaption>Avg. query time over |V(Q)|</figcaption>
</figure>


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/Performance%20(2).png" width="500" alt="Std. dev. of query time over |V(Q)|" />
    <figcaption style="font-style: italic;">Std. dev. of query time over |V(Q)|</figcaption>
</figure>

## Performance of LPI, LPF, and BCPRefine Methods against CECI in the Filtering Step
l2Match derives its filtering strategies from CECI algorithm. Hence, the effectiveness of LPI, LPF, and BCPRefine methods are compared to CECI's filtering performance in terms of the candidate count $\Sigma_{u\in V(Q)}|C(u)|$ and the time taken for the filtering step (filtering time) as shown in figure below. On average, l2Match has 14.39\% lesser candidate count and 39.85\% shorter filtering time than CECI. Average performance different is calculated as the percentage decrease $\overline{p}=100(\Sigma_{i=1}^{N} (x_i^{CECI}-x_i^{l2Match})\div x_i^{CECI})\div N$ where $x$ is the value of candidate count or filtering time, and $N$ is the number of easy (E) queries.


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/l2Match_CECI.png" width="500" alt="Avg. and standard deviation of difference in filtering time and candidate count over the number of query edges |E(Q)|" />
    <figcaption style="font-style: italic;">Avg. and standard deviation of difference in filtering time and candidate count over the number of query edges |E(Q)|</figcaption>
</figure>

## Performance of LPI, and LPF Methods against CECI, CFL, GQL, and DPiso in the Filtering Step
Figures attached below depicts a positive improvement on DPiso in terms of the candidate count, and two obvious positive improvements on GQL and DPiso in terms of the filtering time. The candidate count of the DPiso algorithm, which employs only Label Degree Filtering constraint, is decreased by 16.19\% on average. However, the candidate count of CECI, CFL, and GQL algorithms are barely affected. This proves that the LPF and NLF method have equivalent pruning power. On the other hands, LPI and LPF methods reduces the filtering time by 9.60\% on CFL, 20.75\% on GQL, and 16.94\% on DPiso. Yet, the filtering time of CECI is increased by 0.05\%. We suspect that the overhead to access small amount of neighbors of any vertex in the LPI is greater than accessing the set of all neighbors directly when the number of neighbors $|N(u)|$ of that vertex is relatively small.


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/compareLPF_filter.png" width="500" alt="Avg. and standard deviation of difference in filtering time over the number of query edges |E(Q)|" />
    <figcaption style="font-style: italic;">Avg. and standard deviation of difference in filtering time over the number of query edges |E(Q)|</figcaption>
</figure>


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/compareLPF_cand.png" width="500" alt="Avg. and standard deviation of difference in candidate count over the number of query edges |E(Q)|" />
    <figcaption style="font-style: italic;">Avg. and standard deviation of difference in candidate count over the number of query edges |E(Q)|</figcaption>
</figure>


## Performance of JR Method against CECI, CFL, and GQL in the Enumeration Step
Figures below show the average improvement achieved with JR method in the enumeration step are 12.35\% on CECI, 5.88\% on CFL and 41.62\% on GQL in terms of the query time (total time elapsed), and 46.47\% on CECI, 55.26\% on CFL and 52.94\% on GQL in terms of the number of search nodes. 
Average performance different calculated as the percentage decrease $\overline{p}=100(\Sigma_{i=1}^{N} (x_i^{ORI}-x_i^{JR})\div x_i^{ORI})\div N$ where $x$ is the value of the query time or the number of search nodes, $ORI$ and $JR$ represents the original and optimized algorithm respectively, and $N$ is the number of easy (E) queries.
DPiso is excluded as JR method is not compatible with dynamic enumeration order of the DPiso algorithm. 


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/compareJR_query.png" width="500" alt="Avg. and standard deviation of difference in query time over the number of query edges |E(Q)|" />
    <figcaption style="font-style: italic;">Avg. and standard deviation of difference in query time over the number of query edges |E(Q)|</figcaption>
</figure>


<figure>
    <img src="https://github.com/AltheimIcarus/l2Match/blob/master/Result/compareJR_callcount.png" width="500" alt="Avg. and standard deviation of difference in the number of search nodes over the number of query edges |E(Q)|" />
    <figcaption style="font-style: italic;">Avg. and standard deviation of difference in the number of search nodes over the number of query edges |E(Q)|</figcaption>
</figure>
