# BTechProject
A signed network is a graph containing weighted edges mapped to real numbers, which in turn signify the type and extent of relationships between vertices. The sign denotes the type, whether cohesive or opposite, and the magnitude represents the extent of this relationship. Signed graphs find vital importance as various real life practical scenarios are modelled as signed networks. This project revolves around developing efficient computational methods for finding k sub-graphs in a signed network, such that they are cohesive within themselves, and oppositive among each other. We denote this as k-Opposite Cohesive Groups (k-OCG). The algorithm returns a set of k-subgraphs where edges within each subgraph are dense and cohesive, and these subgraphs have dense opposite edges between them. Thus, they can be seen as groups of like minded vertices against each other. In a large graph, these subgraphs are usually very small, thus finding them becomes a challenge. Many existing algorithms are capable of doing this, but most of them are inefficient to large yet sparse graphs. The FOCG algorithm [https://github.com/lingyangchu/KOCG.SIGKDD2016] outclasses all it's opponents when large and sparse graphs come under consideration, also being the most efficient. But the FOCG algorithm also falls short in spotting 'k' significant OCGs, as it penalizes overlaps between subgraphs without considering their individual contribution to the overall result.
Our algorithm introduces better way for implementation of certain crucial steps in the existing FOCG algorithm, which takes into account the weighted sum of all the subgraphs. This ensures that if no more than k' significant Opposite Cohesive Groups exist, then it returns the most relevant k' groups only on a query for k > k'.
