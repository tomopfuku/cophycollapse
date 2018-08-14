package cophycollapse

import (
	"fmt"
	"math/rand"
)

/*
type HCSearch struct {
	Tree            *Node
	PreorderNodes   []*Node
	Clusters        map[int]*Cluster
	SiteAssignments map[int]int
	Gen             int
	Threads         int
	Workers         int
	ClustOutFile    string
	LogOutFile      string
	K               int
	PrintFreq       int
}
*/
func (s *HCSearch) RunEM() {
	for i := 0; i < s.Gen; i++ {
		if i%s.PrintFreq == 0 {
			fmt.Println("ITERATION", i)
		}
		s.updateClusters()
		s.updateMixtureBranchLengths()
	}
	fmt.Println(len(s.Clusters), s.ClusterString())
}

func (s *HCSearch) SplitEM() {
	clen := len(s.Clusters)
	if clen > 1 {
		for i := 0; i < s.SplitGen; i++ {
			s.updateClusters()
			s.updateMixtureBranchLengths()
			//fmt.Println(i)
			//fmt.Println(s.ClusterString())
		}
		//fmt.Println(clen, s.ClusterString())
	}
}

func (s *HCSearch) updateMixtureBranchLengths() {
	for _, v := range s.Clusters {
		if len(v.Sites) == 0 {
			continue
		}
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.LEN = v.BranchLengths[i]
		}
		//ClusterMissingTraitsEM(s.Tree, v, 100)
		IterateLengthsWeighted(s.Tree, v, 40)
	}
}

func (s *HCSearch) updateClusters() {
	var weights map[int]float64
	for k, v := range s.SiteAssignments {
		weights = s.siteClusterUpdate(k, v)
		for l, c := range s.Clusters {
			c.SiteWeights[k] = weights[l]
			//fmt.Println(k, l, weights[l])
		}
	}
}

func (s *HCSearch) clusterLL() {
	for _, c := range s.Clusters {
		curll := 0.0
		for site := range c.Sites {
			curll += SingleSiteLL(s.Tree, site)
		}
		c.LogLike = curll
	}
}

/*/ClusterString will return a string of the current set of clusters
func (s *HCSearch) ClusterString() string {
	var buffer bytes.Buffer
	cSet := s.Clusters
	for _, like := range cSet {
		buffer.WriteString("(")
		for ind, site := range like.Sites {
			cur := strconv.Itoa(site)
			buffer.WriteString(cur)
			stop := len(like.Sites) - 1
			if ind != stop {
				buffer.WriteString(",")
			}
		}
		buffer.WriteString(");")
	}
	return buffer.String()
}
*/

func (s *HCSearch) siteClusterUpdate(site int, siteClusterLab int) (weights map[int]float64) {
	siteCluster := s.Clusters[siteClusterLab]
	bestLL := -1000000000000.
	var bestClustLab int
	var bestClust *Cluster
	llsum := 0.0
	llmap := make(map[int]float64)
	weights = map[int]float64{}
	for k, v := range s.Clusters {
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.LEN = v.BranchLengths[i]
		}
		curll := SingleSiteLL(s.Tree, site)
		llmap[k] = curll
		llsum += curll

		if curll > bestLL {
			bestLL = curll
			bestClustLab = k
			bestClust = v
		}
	}
	for k, v := range llmap {
		weights[k] = v / llsum
	}

	if bestClustLab != siteClusterLab { //move the site to its new cluster if the best cluster has changed
		var newSlice []int
		for _, s := range siteCluster.Sites {
			if s != site {
				newSlice = append(newSlice, s)
			}
		}
		siteCluster.Sites = newSlice
		bestClust.Sites = append(bestClust.Sites, site)
		s.SiteAssignments[site] = bestClustLab
	}
	return
}

func InitEMSearch(tree *Node, gen int, k int, pr int) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.PreorderNodes = tree.PreorderArray()
	s.Gen = gen
	s.K = k
	//s.startingClustersEMOnly()
	s.singleStartingCluster()
	s.perturbAndUpdate(3)
	s.PrintFreq = pr
	return s
}

/*
this gives starting clusters when K is unknown (for basically a prior-free DPP-style mixture model)
func (search *HCSearch) startingClusters() {
	clus := make(map[int]*Cluster)
	lab := 0
	siteClust := make(map[int]int)
	for k := range search.Tree.CONTRT {
		cur := new(Cluster)
		cur.Sites = append(cur.Sites, k)
		ClusterMissingTraitsEM(search.Tree, cur, 10)
		clus[lab] = cur
		siteClust[k] = lab
		lab++
	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
}
*/

func (search *HCSearch) startingClustersEMOnly() {
	clus := make(map[int]*Cluster)
	siteClust := make(map[int]int)
	var clustLabs []int
	for i := 0; i < search.K; i++ { //create clusters
		cur := new(Cluster)
		clus[i] = cur
		clustLabs = append(clustLabs, i)
		cur.SiteWeights = map[int]float64{}
	}
	for k := range search.Tree.CONTRT {
		lab := rand.Intn(search.K)
		/*
			var lab int
			if k < 50 {
				lab = 0
			} else {
				lab = 1
			}
		*/
		cur := clus[lab]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = lab
	}
	stWt := 1.0 / float64(search.K)
	search.Clusters = clus
	search.SiteAssignments = siteClust
	for _, cur := range search.Clusters {
		if len(cur.Sites) == 0 {
			for range search.PreorderNodes {
				cur.BranchLengths = append(cur.BranchLengths, rand.Float64())
			}
			continue
		}

		for i := range search.SiteAssignments {
			cur.SiteWeights[i] = stWt
		}
		//ClusterMissingTraitsEM(search.Tree, cur, 10)
		IterateLengthsWeighted(search.Tree, cur, 40)
	}

}
