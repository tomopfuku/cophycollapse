package cophycollapse

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"math"
	"math/rand"
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
	TreeOutFile     string
	LogOutFile      string
	K               int
	PrintFreq       int
	CurrentAIC      float64
	NumTraits       float64
	Criterion       int
}

func (s *HCSearch) Run() {
	f, err := os.Create(s.TreeOutFile)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)

	count := 0
	for {
		clcount := len(s.Clusters)
		fmt.Println(clcount)
		if clcount == 1 {
			break
		}
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
	for _, c := range s.Clusters {
		for i, n := range s.PreorderNodes {
			n.LEN = c.BranchLengths[i]
		}
		fmt.Fprint(w, s.Tree.Newick(true)+";"+"\n")
	}
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

func (s *HCSearch) PerturbedRun() {
	f, err := os.Create(s.TreeOutFile)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)

	count := 0
	var bestTrees, bestClust string
	bestAIC := 1000000000.
	for {
		clcount := len(s.Clusters)
		fmt.Println(clcount)

		quit := s.bestClusterJoin()
		if clcount == 1 {
			quit = true
		}
		if quit == true {
			//if math.Abs(bestAIC-s.CurrentAIC) < 4. || count > 1000 {
			if s.CurrentAIC < bestAIC {
				bestAIC = s.CurrentAIC
				bestClust = s.ClusterString()
				bestTrees = ""
				for _, c := range s.Clusters {
					for i, n := range s.PreorderNodes {
						n.LEN = c.BranchLengths[i]
					}
					bestTrees += s.Tree.Newick(true) + ";" + "\n"
				}
			}
			if count > s.Gen {
				break
			}
			fmt.Println(bestAIC, bestClust)
			fmt.Println("Hill climb got stuck. Perturbing the state and trying again to reduce.")
			s.perturbClusters()

		}
		count++
		if count > 10000000000 {
			fmt.Println("couldn't find clusters to join")
			os.Exit(0)
		}
	}
	fmt.Println(bestAIC, bestClust)
	fmt.Fprint(w, bestTrees)
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()

}

func (s *HCSearch) perturbClusters() {
	for i := 0; i < 3; i++ {
		for k, v := range s.SiteAssignments {
			s.siteClusterUpdate(k, v)
		}
		for _, v := range s.Clusters {
			ClusterMissingTraitsEM(s.Tree, v, 10)
		}
	}
	for _, c := range s.Clusters {
		if len(c.Sites) == 0 {
			continue
		}
		clustll := 0.0
		for _, site := range c.Sites {
			cur := SingleSiteLL(s.Tree, site)
			clustll += cur
		}
		c.LogLike = clustll
	}
	s.CurrentAIC = s.calcAIC()
}

func (s *HCSearch) siteClusterUpdate(site int, siteClusterLab int) {
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
	if len(siteCluster.Sites) != 1 {
		var sendsites []int
		sendsites = append(sendsites, site)
		GreedyIterateLengthsMissing(s.Tree, sendsites, 50)
		selfLL := SingleSiteLL(s.Tree, site)
		if selfLL > bestLL {
			bestLL = selfLL
			newLab := MaxClustLab(s.Clusters) + 1
			bestClustLab = newLab
			selfClust := new(Cluster)
			selfClust.LogLike = selfLL //valid because this is a single site cluster
			for _, n := range s.PreorderNodes {
				selfClust.BranchLengths = append(bestClust.BranchLengths, n.LEN)
			}

			s.Clusters[newLab] = selfClust
			//fmt.Println(s.Clusters[newLab])
			bestClust = selfClust
		}
	}
	//fmt.Println(siteClusterLab, bestClustLab)
	if bestClustLab != siteClusterLab { //move the site to its new cluster if the best cluster has changed
		if len(siteCluster.Sites) == 1 { // delete the cluster containing the current site if it is a singleton
			delete(s.Clusters, siteClusterLab)
		} else { // delete the current site from its existing cluster if it shares the cluster with other sites
			var newSlice []int
			for _, s := range siteCluster.Sites {
				if s != site {
					newSlice = append(newSlice, s)
				}
			}
			siteCluster.Sites = newSlice
		}
		//fmt.Println("BEST", bestClust)
		bestClust.Sites = append(bestClust.Sites, site)
		s.SiteAssignments[site] = bestClustLab
	}
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
	var newK int
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
			for _, n := range s.PreorderNodes[1:] {
				r := rand.Float64()
				n.LEN = r
			}
			clen := len(s.Clusters)
			if clen <= 15 {
				GreedyIterateLengthsMissing(s.Tree, proposedSites, 100)
			} else if clen <= 25 && clen > 15 {
				GreedyIterateLengthsMissing(s.Tree, proposedSites, 60)
			} else if clen > 25 {
				GreedyIterateLengthsMissing(s.Tree, proposedSites, 50)
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
				newK = i
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
			s.SiteAssignments[site] = newK
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

func InitGreedyHC(tree *Node, gen int, pr int, crit int, rstart bool, k int, treefl string) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.TreeOutFile = treefl
	s.PreorderNodes = tree.PreorderArray()
	s.Gen = gen
	s.Criterion = crit
	s.K = k
	if rstart == false {
		s.startingClusters()
	} else if rstart == true {
		s.randomStartingClusters()
	}
	s.PrintFreq = pr
	return s
}

func TransferGreedyHC(tree *Node, gen int, pr int, crit int, clus map[int]*Cluster, siteAssign map[int]int, treefl string) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.TreeOutFile = treefl
	s.PreorderNodes = tree.PreorderArray()
	s.Gen = gen
	s.Criterion = crit
	s.PrintFreq = pr
	s.Clusters = clus
	for k, v := range s.Clusters {
		if len(v.Sites) == 0 {
			delete(s.Clusters, k) //delete empty clusters
		}
	}
	s.SiteAssignments = siteAssign
	s.NumTraits = math.Log(float64(len(clus))) //* float64(tipcount)
	s.CurrentAIC = s.calcAIC()
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

func (search *HCSearch) combineAndCalcAIC() {
	newCluster := new(Cluster)
	for k, c := range search.Clusters {
		for _, site := range c.Sites {
			newCluster.Sites = append(newCluster.Sites, site)
		}
		delete(search.Clusters, k)
	}
	MissingTraitsEM(search.Tree, 100)
	ll := CalcUnrootedLogLike(search.Tree, true)
	newCluster.LogLike = ll
	fmt.Println("Single cluster AIC:", search.calcAIC())
}

func (search *HCSearch) randomStartingClusters() {
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
	for k, cur := range clus {
		if len(cur.Sites) != 0 {
			ClusterMissingTraitsEM(search.Tree, cur, 10)
			clustll := 0.0
			for _, site := range cur.Sites {
				curll := SingleSiteLL(search.Tree, site)
				clustll += curll
			}
			cur.LogLike = clustll
		} else {
			delete(clus, k) //delete empty clusters
		}
	}

	search.Clusters = clus
	search.SiteAssignments = siteClust
	search.NumTraits = math.Log(float64(len(clus))) //* float64(tipcount)
	search.CurrentAIC = search.calcAIC()
}
