package cophycollapse

type Cluster struct {
	Sites         []int // this stores all of the sites with a preference for this cluster
	BranchLengths []float64
	LogLike       float64
	SiteWeights   map[int]float64 // this will store the probability that each site in the MATRIX belongs here.
}

func (c *Cluster) CalcLL(tree *Node) {
	c.LogLike = SubUnrootedLogLikeParallel(tree, c.Sites, 6)
}
