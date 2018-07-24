package cophycollapse

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"strconv"
)

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
	CurrentAIC      float64
	NumTraits       float64
	Criterion       int
}

func (s *HCSearch) Run() {
	count := 0
	for {
		fmt.Println(len(s.Clusters))
		quit := s.bestClusterJoin()
		if quit == true {
			break
		}
		count++
		if count > 1000000 {
			fmt.Println("couldn't find clusters to join")
			os.Exit(0)
		}
	}
	fmt.Println(s.ClusterString())
}

func (s *HCSearch) calcAIC() (aic float64) {
	params := 0.
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	ll := 0.
	for _, c := range s.Clusters {
		params += blCount
		ll += c.LogLike
	}
	if s.Criterion == 0 {
		aic = (2. * params) - (2. * ll)
	} else if s.Criterion == 1 {
		aic = (s.NumTraits * params) - (2. * ll)
	}
	//aic += (((2. * (params * params)) + (2. * params)) / (s.NumTraits - params - 1.))

	return
}

func (s *HCSearch) unjoinedLikeParams(exclude1, exclude2 int) (params, ll float64) {
	params = 0.
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	ll = 0.
	for i, c := range s.Clusters {
		if i == exclude1 || i == exclude2 {
			continue
		}
		params += blCount
		ll += c.LogLike
	}
	//aic = (2. * params) - (2. * ll)
	return
}

func (s *HCSearch) bestClusterJoin() (quit bool) {
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	bestAIC := 1000000000000.
	var best1, best2 *Cluster
	var bestll float64
	var bestBranchLengths []float64
	var deleteK int
	for i, c1 := range s.Clusters {
		for j, c2 := range s.Clusters {
			if i >= j {
				continue
			}
			var proposedSites []int
			for _, site := range c1.Sites {
				proposedSites = append(proposedSites, site)
			}
			for _, site := range c2.Sites {
				proposedSites = append(proposedSites, site)
			}
			params, ll := s.unjoinedLikeParams(i, j)
			//fmt.Println(params, len(s.Clusters), len(s.PreorderNodes)-1)
			clen := len(s.Clusters)
			if clen < 11 {
				GreedyIterateLengthsMissing(s.Tree, proposedSites, 50)
			} else if clen < 21 && clen > 11 {
				GreedyIterateLengthsMissing(s.Tree, proposedSites, 20)
			} else if clen > 21 {
				GreedyIterateLengthsMissing(s.Tree, proposedSites, 40)
			}
			clustll := 0.0
			for _, site := range proposedSites {
				cur := SingleSiteLL(s.Tree, site)
				clustll += cur
			}
			ll += clustll
			params += blCount
			aic := 0.
			if s.Criterion == 0 {
				aic = (2. * params) - (2. * ll)
			} else if s.Criterion == 1 {
				aic = (s.NumTraits * params) - (2. * ll)
				//fmt.Println(aic)
			}
			//aic += (((2. * (params * params)) + (2. * params)) / (s.NumTraits - params - 1.))
			if aic < bestAIC {
				//fmt.Println(bestAIC)
				//fmt.Println((((2. * (params * params)) + (2. * params)) / (s.NumTraits - params - 1.)))
				bestAIC = aic
				best1 = c1
				best2 = c2
				deleteK = j
				bestll = clustll
				var bl []float64
				for _, n := range s.PreorderNodes {
					bl = append(bl, n.LEN)
				}
				bestBranchLengths = bl
			}
		}
	}
	fmt.Println(s.CurrentAIC, bestAIC)
	if bestAIC < s.CurrentAIC+10. {
		for _, site := range best2.Sites {
			best1.Sites = append(best1.Sites, site)
		}
		best1.LogLike = bestll
		best1.BranchLengths = bestBranchLengths

		s.CurrentAIC = bestAIC

		delete(s.Clusters, deleteK)
		quit = false
	} else {
		quit = true
	}
	return
}

func InitGreedyHC(tree *Node, gen int, pr int, crit int) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.PreorderNodes = tree.PreorderArray()
	s.Gen = gen
	s.Criterion = crit
	s.startingClusters()
	s.PrintFreq = pr
	return s
}

//this gives starting clusters when K is unknown (for basically a prior-free DPP-style mixture model)
func (search *HCSearch) startingClusters() {
	clus := make(map[int]*Cluster)
	lab := 0
	siteClust := make(map[int]int)
	for k := range search.Tree.CONTRT {
		cur := new(Cluster)
		cur.Sites = append(cur.Sites, k)
		ClusterMissingTraitsEM(search.Tree, cur, 10)
		cur.LogLike = SingleSiteLL(search.Tree, k)
		clus[lab] = cur
		siteClust[k] = lab
		lab++

	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
	search.NumTraits = math.Log(float64(len(clus))) //* float64(tipcount)
	search.CurrentAIC = search.calcAIC()
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.CHLD) == 0 {
			tipcount++
		}
	}
}

//ClusterString will return a string of the current set of clusters
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
