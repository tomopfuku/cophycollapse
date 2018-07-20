package cophycollapse

import (
	"bytes"
	"fmt"
	"math/rand"
	"strconv"
)

type ClustSearch struct {
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

type Cluster struct {
	Sites         []int
	BranchLengths []float64
}

func (s *ClustSearch) Run() {
	for i := 0; i < s.Gen; i++ {
		if i%s.PrintFreq == 0 {
			fmt.Println("ITERATION", i)
		}
		s.updateClusters()
		s.updateClusterBranchLengths()
	}
	fmt.Println(len(s.Clusters), s.ClusterString())
}

func (s *ClustSearch) updateClusterBranchLengths() {
	for _, v := range s.Clusters {
		ClusterMissingTraitsEM(s.Tree, v, 10)
	}
}

func (s *ClustSearch) updateClusters() {
	for k, v := range s.SiteAssignments {
		s.siteClusterUpdate(k, v)
	}
}

//ClusterString will return a string of the current set of clusters
func (s *ClustSearch) ClusterString() string {
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

func (s *ClustSearch) siteClusterUpdate(site int, siteClusterLab int) {
	siteCluster := s.Clusters[siteClusterLab]
	bestLL := -1000000000000.
	var bestClustLab int
	var bestClust *Cluster
	for k, v := range s.Clusters {
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.LEN = v.BranchLengths[i]
		}
		curll := SingleSiteLL(s.Tree, site)
		if curll > bestLL {
			bestLL = curll
			bestClustLab = k
			bestClust = v
		}
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
}

func InitGreedySearch(tree *Node, gen int, k int, pr int) *ClustSearch {
	s := new(ClustSearch)
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
func (search *ClustSearch) startingClusters() {
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

func (search *ClustSearch) startingClusters() {
	clus := make(map[int]*Cluster)
	siteClust := make(map[int]int)
	var clustLabs []int
	for i := 0; i < search.K; i++ { //create clusters
		cur := new(Cluster)
		clus[i] = cur
		clustLabs = append(clustLabs, i)
	}
	for k := range search.Tree.CONTRT {
		lab := rand.Intn(search.K)
		cur := clus[lab]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = lab
	}
	for _, cur := range clus {
		ClusterMissingTraitsEM(search.Tree, cur, 10)
	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
}
