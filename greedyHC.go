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
	RunName         string
	LogOutFile      string
	K               int
	PrintFreq       int
	CurrentAIC      float64
	NumTraits       float64
	Criterion       int
	SavedConfig     []*SiteConfiguration
	CurBestAIC      float64
	JoinLikes       map[int]map[int]float64
	SplitGen        int
}

func (s *HCSearch) NewSiteConfig() *SiteConfiguration {
	config := new(SiteConfiguration)
	var sitemap = map[int]map[int]bool{}
	var treemap = map[int]string{}
	count := 0
	for _, v := range s.Clusters {
		sitemap[count] = map[int]bool{}
		for _, site := range v.Sites {
			sitemap[count][site] = true
		}
		for i, node := range s.PreorderNodes {
			node.LEN = v.BranchLengths[i]
		}
		treemap[count] = s.Tree.Newick(true)
		count++
	}
	config.AIC = s.CurrentAIC
	config.Sites = sitemap
	config.ClusterTrees = treemap
	config.ClusterString = s.ClusterString()
	return config
}

func (s *HCSearch) Run() {
	f, err := os.Create(s.RunName)
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

func (s *HCSearch) CalcRelLikes() (denom float64) {
	denom = 0.
	var deltaAICi float64
	for _, c := range s.SavedConfig {
		if c.AIC == s.CurBestAIC {
			denom += 1.0
			continue
		}
		deltaAICi = c.AIC - s.CurBestAIC
		denom += math.Exp(-0.5 * deltaAICi)
	}
	return
}

func (s *HCSearch) CheckCluster(checkConfig *SiteConfiguration) (keep bool) {
	keep = true
	denom := s.CalcRelLikes()
	relLike := math.Exp(-0.5*checkConfig.AIC - s.CurBestAIC)
	denom += relLike
	weight := relLike / denom
	//fmt.Println("AIC weight:", weight, relLike)
	if weight < .01 {
		keep = false
		return
	}
	for _, c := range s.SavedConfig {
		seen := checkConfig.Equals(c)
		if seen == true {
			keep = false
		}
		return
	}
	return
}

func (s *HCSearch) PerturbedRun() {
	count := 0
	var bestClust string
	bestAIC := 1000000000.
	var bestK int
	fmt.Println("iteration\tcurrent score\tbest score\tbest K")
	for {
		clcount := len(s.Clusters)
		//fmt.Println(clcount)
		var quit bool
		if clcount == 1 {
			quit = true
		} else {
			quit = s.bestClusterJoin()
		}
		if quit == true {
			count++
			config := s.NewSiteConfig()
			var keep bool
			if s.CurrentAIC < bestAIC {
				bestAIC = s.CurrentAIC
				bestK = clcount
				s.CurBestAIC = bestAIC
				bestClust = s.ClusterString()
				keep = true
			} else {
				keep = s.CheckCluster(config)
			}
			if keep == true {
				s.SavedConfig = append(s.SavedConfig, config)
			}
			if count > s.Gen {
				break
			}
			//fmt.Println(bestAIC, bestClust)
			//fmt.Println("Hill climb got stuck. Perturbing the state and trying again to reduce.")
			fmt.Println(count, s.CurrentAIC, bestAIC, bestK)
			s.perturbClusters()
			s.JoinLikes = make(map[int]map[int]float64)
		}

		if count > 10000000000 {
			fmt.Println("couldn't find clusters to join")
			os.Exit(0)
		}
	}
	s.RefineSavedClusterings()
	//fmt.Println(len(s.SavedConfig))
	s.WriteBestClusters()
	s.WriteClusterTrees()
	fmt.Println(bestAIC, bestClust)

}

func (s *HCSearch) WriteBestClusters() {
	f, err := os.Create(s.RunName + "_bestClusters")
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)
	for _, c := range s.SavedConfig {
		fmt.Fprint(w, strconv.FormatFloat(c.AIC, 'f', 6, 64)+"\t"+c.ClusterString+"\n")
	}
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

func (s *HCSearch) WriteClusterTrees() {
	for i, c := range s.SavedConfig {
		stringCount := strconv.Itoa(i)
		f, err := os.Create(s.RunName + "_config" + stringCount + "TREES")
		if err != nil {
			log.Fatal(err)
		}
		w := bufio.NewWriter(f)
		for lab, sitels := range c.Sites {
			var buffer bytes.Buffer
			buffer.WriteString("| ")
			//fmt.Println(sitels)
			for site := range sitels {
				cur := strconv.Itoa(site)
				buffer.WriteString(cur)
				buffer.WriteString(" ")
			}
			buffer.WriteString("| \n\n")
			fmt.Fprint(w, buffer.String())
			fmt.Fprint(w, c.ClusterTrees[lab])
			fmt.Fprint(w, ";\n\n")
		}

		err = w.Flush()
		if err != nil {
			log.Fatal(err)
		}
		f.Close()
	}
}

func (s *HCSearch) RefineSavedClusterings() {
	var kept []*SiteConfiguration
	for _, c := range s.SavedConfig {
		if c.AIC == s.CurBestAIC {
			kept = append(kept, c)
			continue
		}
		denom := s.CalcRelLikes()
		relLike := math.Exp(-0.5*c.AIC - s.CurBestAIC)
		weight := relLike / denom
		if weight > 0.01 {
			kept = append(kept, c)
		}
	}
	s.SavedConfig = kept
}

func (s *HCSearch) perturbClusters() {
	for i := 0; i < 10; i++ {
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
	nclust := len(s.Clusters)
	//llsum := 0.0
	//llmap := make(map[int]float64)
	for k, v := range s.Clusters {
		for i, n := range s.PreorderNodes { //assign current cluster's branch lengths
			n.LEN = v.BranchLengths[i]
		}
		curll := SingleSiteLL(s.Tree, site)
		/*
			llmap[k] = curll
			llsum+=curll
		*/
		if curll > bestLL {
			bestLL = curll
			bestClustLab = k
			bestClust = v
		}
	}
	/*
		for k,v := range llmap {
			fmt.Println(k,v/llsum)
		}
		os.Exit(0)
	*/
	if len(siteCluster.Sites) != 1 && nclust < s.K {
		var sendsites []int
		sendsites = append(sendsites, site)
		for _, n := range s.PreorderNodes[1:] {
			r := rand.Float64()
			n.LEN = r
		}
		GreedyIterateLengthsMissing(s.Tree, sendsites, 10)
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
	} else if s.Criterion == 2 {
		aic = (2. * params) - (2. * ll)
		aic -= ((2. * params) * (params + 1.)) / (s.NumTraits - params - 2.)
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
	//fmt.Println(ll)
	//aic = (2. * params) - (2. * ll)
	return
}

func (s *HCSearch) bestClusterJoin() (quit bool) {
	blCount := float64(len(s.PreorderNodes)) - 1. //subtract 1 because you don't estimate the root node length
	bestAIC := 1000000000000.
	//var best1, best2 *Cluster
	var bestll float64
	var bestBranchLengths []float64
	var bestClusterSites []int
	var deleteK int
	var newK int
	clen := len(s.Clusters)
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
			precalc := true
			var clustll float64
			if _, ok := s.JoinLikes[i][j]; !ok {
				if _, ok := s.JoinLikes[i]; !ok {
					s.JoinLikes[i] = map[int]float64{}
				}
				if _, ok := s.JoinLikes[j][i]; !ok {
					precalc = false // pairwise calc needs to be done
				} else {
					clustll = s.JoinLikes[j][i]
				}
			} else {
				clustll = s.JoinLikes[i][j]
			}
			if precalc == false {
				if clen <= 15 {
					GreedyIterateLengthsMissing(s.Tree, proposedSites, 100)
				} else if clen <= 25 && clen > 15 {
					GreedyIterateLengthsMissing(s.Tree, proposedSites, 60)
				} else if clen > 25 {
					GreedyIterateLengthsMissing(s.Tree, proposedSites, 30)
				}
				clustll = 0.0
				for _, site := range proposedSites {
					cur := SingleSiteLL(s.Tree, site)
					clustll += cur
				}
				s.JoinLikes[i][j] = clustll
			}
			ll += clustll
			params += blCount
			aic := 0.
			if s.Criterion == 0 {
				aic = (2. * params) - (2. * ll)
			} else if s.Criterion == 1 {
				aic = (s.NumTraits * params) - (2. * ll)
				//fmt.Println(aic)
			} else if s.Criterion == 2 {
				aic = (2. * params) - (2. * ll)
				aic -= ((2. * params) * (params + 1.)) / (s.NumTraits - params - 2.)
			}
			if aic < bestAIC {
				//fmt.Println(bestAIC)
				//fmt.Println((((2. * (params * params)) + (2. * params)) / (s.NumTraits - params - 1.)))
				bestAIC = aic
				/*
					best1 = c1
					best2 = c2
				*/
				deleteK = j

				newK = i
				bestll = clustll
				bestClusterSites = proposedSites
				var bl []float64
				for _, n := range s.PreorderNodes {
					bl = append(bl, n.LEN)
				}
				bestBranchLengths = bl
			}
		}
	}
	//fmt.Println(s.CurrentAIC, bestAIC)
	if bestAIC < s.CurrentAIC { //+10. {
		newLab := MaxClustLab(s.Clusters) + 1
		addClust := new(Cluster)
		addClust.Sites = bestClusterSites
		for _, site := range bestClusterSites {
			s.SiteAssignments[site] = newLab
		}
		/*
				for _, site := range best2.Sites {
					best1.Sites = append(best1.Sites, site)
					s.SiteAssignments[site] = newK
				}

			best1.LogLike = bestll
			best1.BranchLengths = bestBranchLengths
		*/
		addClust.LogLike = bestll
		addClust.BranchLengths = bestBranchLengths
		s.CurrentAIC = bestAIC
		delete(s.Clusters, deleteK)
		delete(s.Clusters, newK)
		s.Clusters[newLab] = addClust
		quit = false
	} else {
		quit = true
	}
	return
}

func InitGreedyHC(tree *Node, gen int, pr int, crit int, rstart bool, k int, runName string, splitgen int) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.RunName = runName
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
	s.JoinLikes = make(map[int]map[int]float64)
	s.SplitGen = splitgen
	return s
}

func TransferGreedyHC(tree *Node, gen int, pr int, crit int, clus map[int]*Cluster, siteAssign map[int]int, runName string, splitgen int) *HCSearch {
	s := new(HCSearch)
	s.Tree = tree
	s.RunName = runName
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
	s.JoinLikes = make(map[int]map[int]float64)
	s.SplitGen = splitgen
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
		//ClusterMissingTraitsEM(search.Tree, cur, 10)
		cur.LogLike = SingleSiteLL(search.Tree, k)
		clus[lab] = cur
		siteClust[k] = lab
		lab++

	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.CHLD) == 0 {
			tipcount++
		}
	}
	search.NumTraits = math.Log(float64(len(clus))) * float64(tipcount)
	search.CurrentAIC = search.calcAIC()
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
		buffer.WriteString(")")
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

func (search *HCSearch) singleStartingCluster() {
	clus := make(map[int]*Cluster)
	siteClust := make(map[int]int)
	var clustLabs []int
	cur := new(Cluster)
	clus[0] = cur
	clustLabs = append(clustLabs, 0)
	for k := range search.Tree.CONTRT {
		//cur := clus[0]
		cur.Sites = append(cur.Sites, k)
		siteClust[k] = 0
	}
	if len(cur.Sites) != 0 {
		//ClusterMissingTraitsEM(search.Tree, cur, 100)
		clustll := 0.0
		for _, site := range cur.Sites {
			curll := SingleSiteLL(search.Tree, site)
			clustll += curll
		}
		cur.LogLike = clustll
	}
	search.Clusters = clus
	search.SiteAssignments = siteClust
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.CHLD) == 0 {
			tipcount++
		}
	}
	search.NumTraits = math.Log(float64(len(clus))) * float64(tipcount)
	search.CurrentAIC = search.calcAIC()
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
			//ClusterMissingTraitsEM(search.Tree, cur, 100)
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
	tipcount := 0
	for _, n := range search.PreorderNodes {
		if len(n.CHLD) == 0 {
			tipcount++
		}
	}
	search.NumTraits = math.Log(float64(len(clus))) * float64(tipcount)
	search.CurrentAIC = search.calcAIC()
}
