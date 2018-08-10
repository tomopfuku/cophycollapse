package cophycollapse

import (
	"bytes"
	"fmt"
	"math/rand"
	"strconv"
)

type EMClustSearch struct {
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

func (s *EMClustSearch) Run() {
	for i := 0; i < s.Gen; i++ {
		if i%s.PrintFreq == 0 {
			fmt.Println("ITERATION", i)
		}
		s.updateClusters()
		s.updateClusterBranchLengths()
	}
	fmt.Println(len(s.Clusters), s.ClusterString())
}

func (s *EMClustSearch) updateClusterBranchLengths() {
	for _, v := range s.Clusters {
		if len(v.Sites) == 0 {
			continue
		}
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.LEN = v.BranchLengths[i]
		}
		//ClusterMissingTraitsEM(s.Tree, v, 100)
		IterateLengthsWeighted(s.Tree, v, 100)
	}
}

func (s *EMClustSearch) updateClusters() (weights map[int]float64) {
	for k, v := range s.SiteAssignments {
		weights = s.siteClusterUpdate(k, v)
		for l, c := range s.Clusters {
			c.SiteWeights[k] = weights[l]
			fmt.Println(k, l, weights[l])
		}
	}
	return
}

//ClusterString will return a string of the current set of clusters
func (s *EMClustSearch) ClusterString() string {
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

func (s *EMClustSearch) siteClusterUpdate(site int, siteClusterLab int) (weights map[int]float64) {
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
		//fmt.Println(k, v/llsum)
	}

	//fmt.Println(siteClusterLab, bestClustLab)
	if bestClustLab != siteClusterLab { //move the site to its new cluster if the best cluster has changed
		//if len(siteCluster.Sites) == 1 { // delete the cluster containing the current site if it is a singleton
		//delete(s.Clusters, siteClusterLab)
		//} else { // delete current site from its cluster
		//if len(siteCluster.Sites) != 1 { // delete the current site from its existing cluster if it shares the cluster with other sites
		var newSlice []int
		for _, s := range siteCluster.Sites {
			if s != site {
				newSlice = append(newSlice, s)
			}
		}
		siteCluster.Sites = newSlice
		//}
		//fmt.Println("BEST", bestClust)
		bestClust.Sites = append(bestClust.Sites, site)
		s.SiteAssignments[site] = bestClustLab
	}
	return
}

func InitEMSearch(tree *Node, gen int, k int, pr int) *EMClustSearch {
	s := new(EMClustSearch)
	s.Tree = tree
	s.PreorderNodes = tree.PreorderArray()
	s.Gen = gen
	s.K = k
	s.startingClusters()
	s.PrintFreq = pr
	return s
}

/*
this gives starting clusters when K is unknown (for basically a prior-free DPP-style mixture model)
func (search *EMClustSearch) startingClusters() {
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

func (search *EMClustSearch) startingClusters() {
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
		cur := clus[lab]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = lab
	}
	for _, cur := range clus {
		//if len(cur.Sites) == 0 {
		for range search.PreorderNodes {
			cur.BranchLengths = append(cur.BranchLengths, rand.Float64())
		}
		for i := range search.SiteAssignments {
			cur.SiteWeights[i] = 0.5
		}
		//continue
		//}
		//ClusterMissingTraitsEM(search.Tree, cur, 10)
	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
}
