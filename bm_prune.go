package cophycollapse

import (
	"fmt"
	"math"
	"os"
)

func postorder(curnode *Node) {
	for _, chld := range curnode.CHLD {
		postorder(chld)
	}
	fmt.Println(curnode.NAME, curnode.CONTRT)
}

//BMPruneRooted will prune BM branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRooted(n *Node) {
	for _, chld := range n.CHLD {
		BMPruneRooted(chld)
	}
	n.PRNLEN = n.LEN
	nchld := len(n.CHLD)
	if nchld != 0 { //&& n.MRK == false {
		var tempChar float64
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		bot := ((1.0 / c0.PRNLEN) + (1.0 / c1.PRNLEN))
		n.PRNLEN += 1.0 / bot
		for i := range n.CHLD[0].CONTRT {
			tempChar = (((1 / c0.PRNLEN) * c1.CONTRT[i]) + ((1 / c1.PRNLEN) * c0.CONTRT[i])) / bot
			n.CONTRT[i] = tempChar
		}
	}
}

//BMPruneRootedSub will prune BM branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRootedSub(n *Node, sites []int) {
	for _, chld := range n.CHLD {
		BMPruneRootedSub(chld, sites)
	}
	n.PRNLEN = n.LEN
	nchld := len(n.CHLD)
	if nchld != 0 { //&& n.MRK == false {
		var tempChar float64
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		bot := ((1.0 / c0.PRNLEN) + (1.0 / c1.PRNLEN))
		n.PRNLEN += 1.0 / bot
		if math.IsNaN(n.PRNLEN) {
			fmt.Println(c0.NAME, c0.PRNLEN, c1.NAME, c1.PRNLEN, sites)
			os.Exit(0)
		}
		for _, i := range sites { //n.CHLD[0].CONTRT {
			tempChar = (((1 / c0.PRNLEN) * c1.CONTRT[i]) + ((1 / c1.PRNLEN) * c0.CONTRT[i])) / bot
			n.CONTRT[i] = tempChar
		}
	}
}

//BMPruneRootedSingle will prune BM branch lens and calculate PIC of a single trait down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRootedSingle(n *Node, i int) {
	for _, chld := range n.CHLD {
		BMPruneRootedSingle(chld, i)
	}
	n.PRNLEN = n.LEN
	nchld := len(n.CHLD)
	if nchld != 0 { //&& n.MRK == false {
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		bot := ((1.0 / c0.PRNLEN) + (1.0 / c1.PRNLEN))
		n.PRNLEN += 1.0 / bot
		if n.PRNLEN == 0. {
			n.PRNLEN = 0.0001
		}
		tempCharacter := (((1 / c0.PRNLEN) * c1.CONTRT[i]) + ((1 / c1.PRNLEN) * c0.CONTRT[i])) / bot
		if math.IsNaN(tempCharacter) {
			fmt.Println(c0.PRNLEN, c1.PRNLEN, c0.CONTRT[i], c1.CONTRT[i], c0.NAME, c1.NAME)
			fmt.Println("you're encountering zero branch lengths while pruning to the root on trait", i, "in the matrix")
			os.Exit(0)
		}
		n.CONTRT[i] = tempCharacter
	}
}

func debugParChld(tree *Node) {
	for _, ch := range tree.CHLD {
		fmt.Println(ch.NAME, "\t")
	}
}

//AssertUnrootedTree is a quick check to make sure the tree passed is unrooted
func AssertUnrootedTree(tree *Node) {
	if len(tree.CHLD) != 3 {
		fmt.Print("BRANCH LENGTHS MUST BE ITERATED ON AN UNROOTED TREE. THIS TREE IS ROOTED.")
		os.Exit(0)
	}
}

//IterateBMLengths will iteratively calculate the ML branch lengths for a particular topology, assuming that traits are fully sampled at the tips. should use MissingTraitsEM if there are missing sites.
func IterateBMLengths(tree *Node, niter int) {
	AssertUnrootedTree(tree)
	itercnt := 0
	for {
		calcBMLengths(tree)
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

//MissingTraitsEM will iteratively calculate the ML branch lengths for a particular topology
func MissingTraitsEM(tree *Node, niter int) {
	AssertUnrootedTree(tree)
	//nodes := tree.PreorderArray()
	//InitMissingValues(nodes)
	itercnt := 0
	for {
		CalcExpectedTraits(tree) //calculate Expected trait values
		calcBMLengths(tree)      //maximize likelihood of branch lengths
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

//GreedyIterateLengthsMissing will iteratively calculate the ML branch lengths for a particular topology and cluster when doing the greedy site clustering procedure.
func GreedyIterateLengthsMissing(tree *Node, sites []int, niter int) {
	AssertUnrootedTree(tree)
	//nodes := tree.PreorderArray()
	//InitMissingValues(nodes)
	rnodes := tree.PreorderArray()
	for i := 0; i < niter; i++ {
		CalcExpectedTraitsSub(tree, sites)         //calculate Expected trait values
		calcBMLengthsSubSites(tree, rnodes, sites) //maximize likelihood of branch lengths
	}
}

//IterateLengthsWeighted will iteratively calculate the ML branch lengths for a particular topology and cluster when doing the greedy site clustering procedure.
func IterateLengthsWeighted(tree *Node, cluster *Cluster, niter int) {
	AssertUnrootedTree(tree)
	//nodes := tree.PreorderArray()
	//InitMissingValues(nodes)
	rnodes := tree.PreorderArray()
	for i := 0; i < niter; i++ {
		CalcExpectedTraits(tree)                                 //calculate Expected trait values
		calcBMLengthsWeighted(tree, rnodes, cluster.SiteWeights) //maximize likelihood of branch lengths
	}
	for _, n := range rnodes {
		cluster.BranchLengths = append(cluster.BranchLengths, n.LEN) //store the newly calculated branch lengths
	}
}

//ClusterMissingTraitsEM will iteratively calculate the ML branch lengths for a particular topology and cluster when doing the greedy site clustering procedure.
func ClusterMissingTraitsEM(tree *Node, cluster *Cluster, niter int) {
	AssertUnrootedTree(tree)
	//nodes := tree.PreorderArray()
	//InitMissingValues(nodes)
	rnodes := tree.PreorderArray()
	for i := 0; i < niter; i++ {
		CalcExpectedTraitsSub(tree, cluster.Sites)         //calculate Expected trait values
		calcBMLengthsSubSites(tree, rnodes, cluster.Sites) //maximize likelihood of branch lengths
	}
	for _, n := range rnodes {
		cluster.BranchLengths = append(cluster.BranchLengths, n.LEN) //store the newly calculated branch lengths
	}
}

//calcBMLengthsWeighted will perform a single pass of the branch length ML estimation using only the sites indicated in sites
func calcBMLengthsWeighted(tree *Node, rnodes []*Node, weights map[int]float64) {
	lnode := 0
	for ind, newroot := range rnodes {
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.CHLD {
			BMPruneRooted(cn)
		}
		TritomyWeightedML(tree, weights)
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}

//calcBMLengths will perform a single pass of the branch length ML estimation using only the sites indicated in sites
func calcBMLengthsSubSites(tree *Node, rnodes []*Node, sites []int) {
	lnode := 0
	for ind, newroot := range rnodes {
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.CHLD {
			BMPruneRootedSub(cn, sites)
		}
		TritomySubML(tree, sites)
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}

//calcBMLengths will perform a single pass of the branch length ML estimation
func calcBMLengths(tree *Node) {
	rnodes := tree.PreorderArray()
	lnode := 0
	for ind, newroot := range rnodes {
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.CHLD {
			BMPruneRooted(cn)
		}
		TritomyML(tree)
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}

//PruneToStar will prune brlens and traits to a root
func PruneToStar(tree *Node) {
	for _, cn := range tree.CHLD {
		BMPruneRooted(cn)
	}
}

//CalcUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that no sites have missing data.
func CalcUnrootedLogLike(tree *Node, startFresh bool) (chll float64) {
	chll = 0.0
	for _, ch := range tree.CHLD {
		curlike := 0.0
		CalcRootedLogLike(ch, &curlike, startFresh)
		chll += curlike
	}
	sitelikes := 0.0
	var tmpll float64
	var contrast, curVar float64
	for i := range tree.CHLD[0].CONTRT {
		tmpll = 0.
		if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false { //do the standard calculation when no subtrees have missing traits
			contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
			curVar = tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
			tmpPRNLEN := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].PRNLEN) / (tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN))
			tmpChar := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].CONTRT[i]) + (tree.CHLD[1].PRNLEN * tree.CHLD[0].CONTRT[i])) / curVar
			contrast = tmpChar - tree.CHLD[2].CONTRT[i]
			curVar = tree.CHLD[2].PRNLEN + tmpPRNLEN
			tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		}
		sitelikes += tmpll
	}
	chll += sitelikes
	return
}

//CalcRootedLogLike will return the BM likelihood of a tree assuming that no data are missing from the tips.
func CalcRootedLogLike(n *Node, nlikes *float64, startFresh bool) {
	for _, chld := range n.CHLD {
		CalcRootedLogLike(chld, nlikes, startFresh)
	}
	nchld := len(n.CHLD)
	if n.MRK == true && startFresh == false {
		if nchld != 0 {
			for _, l := range n.LL {
				*nlikes += l
			}
		}
	} else if n.MRK == false || startFresh == true {
		n.PRNLEN = n.LEN
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
			}
			c0 := n.CHLD[0]
			c1 := n.CHLD[1]
			curlike := float64(0.0)
			var tempChar float64
			for i := range n.CHLD[0].CONTRT {
				curVar := c0.PRNLEN + c1.PRNLEN
				contrast := c0.CONTRT[i] - c1.CONTRT[i]
				curlike += ((-0.5) * ((math.Log(2 * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
				n.LL[i] = curlike
				tempChar = ((c0.PRNLEN * c1.CONTRT[i]) + (c1.PRNLEN * c0.CONTRT[i])) / (curVar)
				n.CONTRT[i] = tempChar
			}
			*nlikes += curlike
			tempBranchLength := n.LEN + ((c0.PRNLEN * c1.PRNLEN) / (c0.PRNLEN + c1.PRNLEN)) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.PRNLEN = tempBranchLength                                                     // need to calculate the "prune length" by adding the length to the uncertainty

			n.MRK = true
		}
	}
}

//SitewiseLogLike will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func SitewiseLogLike(tree *Node) (sitelikes []float64) {
	var tmpll float64
	ch1 := tree.CHLD[0]                     //.PostorderArray()
	ch2 := tree.CHLD[1]                     //.PostorderArray()
	ch3 := tree.CHLD[2]                     //.PostorderArray()
	for site := range tree.CHLD[0].CONTRT { //calculate log likelihood at each site
		tmpll = 0.
		calcRootedSiteLL(ch1, &tmpll, true, site)
		calcRootedSiteLL(ch2, &tmpll, true, site)
		calcRootedSiteLL(ch3, &tmpll, true, site)
		tmpll += calcUnrootedSiteLL(tree, site)
		//fmt.Println(site, tmpll)
		sitelikes = append(sitelikes, tmpll)
	}
	return
}

/*/WeightedUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that no sites have missing data.
func WeightedUnrootedLogLike(tree *Node, startFresh bool, weights []float64) (chll float64) {
	chll = 0.0
	for _, ch := range tree.CHLD {
		n := ch.PostorderArray()
		for i := range ch.CONTRT {
			curlike := calcRootedSiteLL(n, startFresh, i)
			chll += curlike
		}
	}
	sitelikes := 0.0
	var tmpll float64
	var contrast, curVar float64
	for i := range tree.CHLD[0].CONTRT {
		tmpll = 0.
		if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false { //do the standard calculation when no subtrees have missing traits
			contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
			curVar = tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
			tmpPRNLEN := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].PRNLEN) / (tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN))
			tmpChar := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].CONTRT[i]) + (tree.CHLD[1].PRNLEN * tree.CHLD[0].CONTRT[i])) / curVar
			contrast = tmpChar - tree.CHLD[2].CONTRT[i]
			curVar = tree.CHLD[2].PRNLEN + tmpPRNLEN
			tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		}
		sitelikes += tmpll * weights[i]
	}
	chll += sitelikes
	return
}
*/

//MarkAll will mark all of the nodes in a tree ARRAY
func MarkAll(nodes []*Node) {
	for _, n := range nodes {
		n.MRK = true
	}
}

//calcUnrootedNodeLikes  will calculate the likelihood of an unrooted tree at each site (i) of the continuous character alignment
func calcUnrootedSiteLLParallel(tree *Node, i int) (tmpll float64) {
	var contrast, curVar float64
	log2pi := 1.8378770664093453
	if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false { //do the standard calculation when no subtrees have missing traits
		contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
		curVar = tree.CHLD[0].CONPRNLEN[i] + tree.CHLD[1].CONPRNLEN[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		tmpCONPRNLEN := ((tree.CHLD[0].CONPRNLEN[i] * tree.CHLD[1].CONPRNLEN[i]) / (tree.CHLD[0].CONPRNLEN[i] + tree.CHLD[1].CONPRNLEN[i]))
		tmpChar := ((tree.CHLD[0].CONPRNLEN[i] * tree.CHLD[1].CONTRT[i]) + (tree.CHLD[1].CONPRNLEN[i] * tree.CHLD[0].CONTRT[i])) / curVar
		contrast = tmpChar - tree.CHLD[2].CONTRT[i]
		curVar = tree.CHLD[2].CONPRNLEN[i] + tmpCONPRNLEN
		tmpll += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == true { // do standard "rooted" calculation on CHLD[0] and CHLD [1] if CHLD[2] is missing
		contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
		curVar = tree.CHLD[0].CONPRNLEN[i] + tree.CHLD[1].CONPRNLEN[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.CHLD[0].MIS[i] == false && tree.CHLD[2].MIS[i] == false && tree.CHLD[1].MIS[i] == true { // do standard "rooted" calculation on CHLD[0] and CHLD [2] if CHLD[1] is missing
		contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[2].CONTRT[i]
		curVar = tree.CHLD[0].CONPRNLEN[i] + tree.CHLD[2].CONPRNLEN[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false && tree.CHLD[0].MIS[i] == true { // do standard "rooted" calculation on CHLD[1] and CHLD [2] if CHLD[0] is missing
		contrast = tree.CHLD[1].CONTRT[i] - tree.CHLD[2].CONTRT[i]
		curVar = tree.CHLD[1].CONPRNLEN[i] + tree.CHLD[2].CONPRNLEN[i]
		tmpll = ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	}
	return
}

//calcRootedSiteLL will return the BM likelihood of a tree assuming that no data are missing from the tips.
func calcRootedSiteLLParallel(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.CHLD {
		calcRootedSiteLLParallel(chld, nlikes, startFresh, site)
	}
	nchld := len(n.CHLD)
	if n.MRK == true {
		if startFresh == false {
			if nchld != 0 {
				*nlikes += n.LL[site]
			}
		}
	}
	if n.MRK == false || startFresh == true {
		n.CONPRNLEN[site] = n.LEN
		log2pi := 1.8378770664093453
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
				os.Exit(0)
			}
			c0 := n.CHLD[0]
			c1 := n.CHLD[1]
			curlike := float64(0.0)
			var tempChar float64
			curVar := c0.CONPRNLEN[site] + c1.CONPRNLEN[site]
			contrast := c0.CONTRT[site] - c1.CONTRT[site]
			curlike += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = ((c0.CONPRNLEN[site] * c1.CONTRT[site]) + (c1.CONPRNLEN[site] * c0.CONTRT[site])) / (curVar)
			n.CONTRT[site] = tempChar
			*nlikes += curlike
			tempBranchLength := n.CONPRNLEN[site] + ((c0.CONPRNLEN[site] * c1.CONPRNLEN[site]) / (c0.CONPRNLEN[site] + c1.CONPRNLEN[site])) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.CONPRNLEN[site] = tempBranchLength                                                                                            // need to calculate the "prune length" by adding the length to the uncertainty
			n.LL[site] = curlike
			//n.MRK = true
		}
	}
}

func siteTreeLikeParallel(tree, ch1, ch2, ch3 *Node, startFresh bool, weights []float64, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		calcRootedSiteLLParallel(ch1, &tmpll, startFresh, site)
		calcRootedSiteLLParallel(ch2, &tmpll, startFresh, site)
		calcRootedSiteLLParallel(ch3, &tmpll, startFresh, site)
		tmpll += calcUnrootedSiteLLParallel(tree, site)
		tmpll = tmpll * weights[site]
		results <- tmpll
	}
}

//WeightedUnrootedLogLikeParallel will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func WeightedUnrootedLogLikeParallel(tree *Node, startFresh bool, weights []float64, workers int) (sitelikes float64) {
	nsites := len(tree.CHLD[0].CONTRT)
	ch1 := tree.CHLD[0] //.PostorderArray()
	ch2 := tree.CHLD[1] //.PostorderArray()
	ch3 := tree.CHLD[2] //.PostorderArray()
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for w := 0; w < workers; w++ {
		go siteTreeLikeParallel(tree, ch1, ch2, ch3, startFresh, weights, jobs, results)
	}

	for site := 0; site < nsites; site++ {
		jobs <- site
	}
	close(jobs)

	for site := 0; site < nsites; site++ {
		sitelikes += <-results
	}
	return
}

//SingleSiteLL will return the likelihood of a single site
func SingleSiteLL(tree *Node, site int) (sitelike float64) {
	sitelike = 0.0
	var tmpll float64
	ch1 := tree.CHLD[0] //.PostorderArray()
	ch2 := tree.CHLD[1] //.PostorderArray()
	ch3 := tree.CHLD[2] //.PostorderArray()
	tmpll = 0.
	calcRootedSiteLL(ch1, &tmpll, true, site)
	calcRootedSiteLL(ch2, &tmpll, true, site)
	calcRootedSiteLL(ch3, &tmpll, true, site)
	tmpll += calcUnrootedSiteLL(tree, site)
	//fmt.Println(tmpll)
	sitelike = tmpll
	return
}

//WeightedUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This can be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func WeightedUnrootedLogLike(tree *Node, startFresh bool, weights []float64) (sitelikes float64) {
	sitelikes = 0.0
	var tmpll float64
	ch1 := tree.CHLD[0]                     //.PostorderArray()
	ch2 := tree.CHLD[1]                     //.PostorderArray()
	ch3 := tree.CHLD[2]                     //.PostorderArray()
	for site := range tree.CHLD[0].CONTRT { //calculate log likelihood at each site
		tmpll = 0.
		calcRootedSiteLL(ch1, &tmpll, startFresh, site)
		calcRootedSiteLL(ch2, &tmpll, startFresh, site)
		calcRootedSiteLL(ch3, &tmpll, startFresh, site)
		tmpll += calcUnrootedSiteLL(tree, site)
		tmpll = tmpll * weights[site]
		//fmt.Println(tmpll)
		sitelikes += tmpll
	}
	return
}

//MissingUnrootedLogLike will calculate the log-likelihood of an unrooted tree. It is deprecated and is only a convienience for testing LL calculation.
func MissingUnrootedLogLike(tree *Node, startFresh bool) (sitelikes float64) {
	sitelikes = 0.0
	var tmpll float64
	ch1 := tree.CHLD[0]                     //.PostorderArray()
	ch2 := tree.CHLD[1]                     //.PostorderArray()
	ch3 := tree.CHLD[2]                     //.PostorderArray()
	for site := range tree.CHLD[0].CONTRT { //calculate log likelihood at each site
		tmpll = 0.
		calcRootedSiteLL(ch1, &tmpll, startFresh, site)
		calcRootedSiteLL(ch2, &tmpll, startFresh, site)
		calcRootedSiteLL(ch3, &tmpll, startFresh, site)
		tmpll += calcUnrootedSiteLL(tree, site)
		sitelikes += tmpll
	}
	return
}

//calcUnrootedNodeLikes  will calculate the likelihood of an unrooted tree at each site (i) of the continuous character alignment
func calcUnrootedSiteLL(tree *Node, i int) (tmpll float64) {
	var contrast, curVar float64
	if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false { //do the standard calculation when no subtrees have missing traits
		contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
		curVar = tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		tmpPRNLEN := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].PRNLEN) / (tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN))
		tmpChar := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].CONTRT[i]) + (tree.CHLD[1].PRNLEN * tree.CHLD[0].CONTRT[i])) / curVar
		contrast = tmpChar - tree.CHLD[2].CONTRT[i]
		curVar = tree.CHLD[2].PRNLEN + tmpPRNLEN
		tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == true { // do standard "rooted" calculation on CHLD[0] and CHLD [1] if CHLD[2] is missing
		contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
		curVar = tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.CHLD[0].MIS[i] == false && tree.CHLD[2].MIS[i] == false && tree.CHLD[1].MIS[i] == true { // do standard "rooted" calculation on CHLD[0] and CHLD [2] if CHLD[1] is missing
		contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[2].CONTRT[i]
		curVar = tree.CHLD[0].PRNLEN + tree.CHLD[2].PRNLEN
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	} else if tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false && tree.CHLD[0].MIS[i] == true { // do standard "rooted" calculation on CHLD[1] and CHLD [2] if CHLD[0] is missing
		contrast = tree.CHLD[1].CONTRT[i] - tree.CHLD[2].CONTRT[i]
		curVar = tree.CHLD[1].PRNLEN + tree.CHLD[2].PRNLEN
		tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
	}
	return
}

//MissingRootedLogLike will return the BM log-likelihood of a tree for a single site, pruning tips that have missing data
func MissingRootedLogLike(n *Node, startFresh bool) (sitelikes float64) {
	//nodes := n.PostorderArray()
	sitelikes = 0.
	for site := range n.CONTRT {
		calcRootedSiteLL(n, &sitelikes, startFresh, site)
	}
	return
}

//calcRootedSiteLL will return the BM likelihood of a tree assuming that no data are missing from the tips.
func calcRootedSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.CHLD {
		calcRootedSiteLL(chld, nlikes, startFresh, site)
	}
	nchld := len(n.CHLD)
	if n.MRK == true {
		if startFresh == false {
			if nchld != 0 {
				*nlikes += n.LL[site]
			}
		}
	}
	if n.MRK == false || startFresh == true {
		n.PRNLEN = n.LEN
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
				os.Exit(0)
			}
			c0 := n.CHLD[0]
			c1 := n.CHLD[1]
			curlike := float64(0.0)
			var tempChar float64
			curVar := c0.PRNLEN + c1.PRNLEN
			contrast := c0.CONTRT[site] - c1.CONTRT[site]
			curlike += ((-0.5) * ((math.Log(2 * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = ((c0.PRNLEN * c1.CONTRT[site]) + (c1.PRNLEN * c0.CONTRT[site])) / (curVar)
			n.CONTRT[site] = tempChar
			*nlikes += curlike
			tempBranchLength := n.LEN + ((c0.PRNLEN * c1.PRNLEN) / (c0.PRNLEN + c1.PRNLEN)) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.PRNLEN = tempBranchLength                                                     // need to calculate the "prune length" by adding the length to the uncertainty
			n.LL[site] = curlike
			//n.MRK = true
		}
	}
}

//TritomySubML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree using only the sites indicated in sites
func TritomySubML(tree *Node, sites []int) {
	ntraits := len(sites)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for _, i := range sites {
		x1 = tree.CHLD[0].CONTRT[i]
		x2 = tree.CHLD[1].CONTRT[i]
		x3 = tree.CHLD[2].CONTRT[i]
		sumV1 += ((x1 - x2) * (x1 - x3))
		sumV2 += ((x2 - x1) * (x2 - x3))
		sumV3 += ((x3 - x1) * (x3 - x2))
	}
	if sumV1 < 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for _, i := range sites {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV2 += (x1 - x2) * (x1 - x2)
			sumV3 += (x1 - x3) * (x1 - x3)
		}
	} else if sumV2 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for _, i := range sites {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x2 - x1) * (x2 - x1)
			sumV3 += (x2 - x3) * (x2 - x3)
		}
	} else if sumV3 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for _, i := range sites {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	sumV1 = sumV1 / fntraits
	sumV2 = sumV2 / fntraits
	sumV3 = sumV3 / fntraits
	sumV1 = sumV1 - (tree.CHLD[0].PRNLEN - tree.CHLD[0].LEN)
	sumV2 = sumV2 - (tree.CHLD[1].PRNLEN - tree.CHLD[1].LEN)
	sumV3 = sumV3 - (tree.CHLD[2].PRNLEN - tree.CHLD[2].LEN)

	if sumV1 <= 0. {
		sumV1 = 0.0001
	}
	if sumV2 <= 0. {
		sumV2 = 0.0001
	}
	if sumV3 <= 0. {
		sumV3 = 0.0001
	}
	tree.CHLD[0].LEN = sumV1
	tree.CHLD[1].LEN = sumV2
	tree.CHLD[2].LEN = sumV3
}

//TritomyWeightedML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func TritomyWeightedML(tree *Node, weights map[int]float64) {
	//ntraits := len(tree.CHLD[0].CONTRT)
	fntraits := 0.
	for _, w := range weights {
		fntraits += w
	}
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	wt := 0.0
	for i := range tree.CHLD[0].CONTRT {
		wt = weights[i]
		if wt == 0.0 {
			continue
		}
		x1 = tree.CHLD[0].CONTRT[i]
		x2 = tree.CHLD[1].CONTRT[i]
		x3 = tree.CHLD[2].CONTRT[i]

		sumV1 += ((x1 - x2) * (x1 - x3)) * wt
		sumV2 += ((x2 - x1) * (x2 - x3)) * wt
		sumV3 += ((x3 - x1) * (x3 - x2)) * wt
	}
	if sumV1 <= 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range tree.CHLD[0].CONTRT {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV2 += (x1 - x2) * (x1 - x2) * wt
			sumV3 += (x1 - x3) * (x1 - x3) * wt
		}
	} else if sumV2 <= 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := range tree.CHLD[0].CONTRT {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}

			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x2 - x1) * (x2 - x1) * wt
			sumV3 += (x2 - x3) * (x2 - x3) * wt
		}
	} else if sumV3 <= 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := range tree.CHLD[0].CONTRT {
			wt = weights[i]
			if wt == 0.0 {
				continue
			}

			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x3 - x1) * (x3 - x1) * wt
			sumV2 += (x3 - x2) * (x3 - x2) * wt
		}
	}
	sumV1 = sumV1 / fntraits
	sumV2 = sumV2 / fntraits
	sumV3 = sumV3 / fntraits
	sumV1 = sumV1 - (tree.CHLD[0].PRNLEN - tree.CHLD[0].LEN)
	sumV2 = sumV2 - (tree.CHLD[1].PRNLEN - tree.CHLD[1].LEN)
	sumV3 = sumV3 - (tree.CHLD[2].PRNLEN - tree.CHLD[2].LEN)
	if sumV1 <= 0. {
		sumV1 = 0.0001
	}
	if sumV2 <= 0. {
		sumV2 = 0.0001
	}
	if sumV3 <= 0. {
		sumV3 = 0.0001
	}
	//fmt.Println(tree.CHLD[0].NAME, sumV1, tree.CHLD[1].NAME, sumV2, tree.CHLD[2].NAME, sumV3)
	if math.IsNaN(sumV1) || math.IsNaN(sumV2) || math.IsNaN(sumV3) {
		fmt.Println(tree.CHLD[0].NAME, sumV1, tree.CHLD[1].NAME, sumV2, tree.CHLD[2].NAME, sumV3)
		os.Exit(0)
	}
	tree.CHLD[0].LEN = sumV1
	tree.CHLD[1].LEN = sumV2
	tree.CHLD[2].LEN = sumV3
}

//TritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func TritomyML(tree *Node) {
	ntraits := len(tree.CHLD[0].CONTRT)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := range tree.CHLD[0].CONTRT {
		x1 = tree.CHLD[0].CONTRT[i]
		x2 = tree.CHLD[1].CONTRT[i]
		x3 = tree.CHLD[2].CONTRT[i]
		sumV1 += ((x1 - x2) * (x1 - x3))
		sumV2 += ((x2 - x1) * (x2 - x3))
		sumV3 += ((x3 - x1) * (x3 - x2))
	}
	if sumV1 < 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV2 += (x1 - x2) * (x1 - x2)
			sumV3 += (x1 - x3) * (x1 - x3)
		}
	} else if sumV2 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x2 - x1) * (x2 - x1)
			sumV3 += (x2 - x3) * (x2 - x3)
		}
	} else if sumV3 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	sumV1 = sumV1 / fntraits
	sumV2 = sumV2 / fntraits
	sumV3 = sumV3 / fntraits
	sumV1 = sumV1 - (tree.CHLD[0].PRNLEN - tree.CHLD[0].LEN)
	sumV2 = sumV2 - (tree.CHLD[1].PRNLEN - tree.CHLD[1].LEN)
	sumV3 = sumV3 - (tree.CHLD[2].PRNLEN - tree.CHLD[2].LEN)
	if sumV1 <= 0. {
		sumV1 = 0.0001
	}
	if sumV2 <= 0. {
		sumV2 = 0.0001
	}
	if sumV3 <= 0. {
		sumV3 = 0.0001
	}
	//fmt.Println(tree.CHLD[0].NAME, sumV1, tree.CHLD[1].NAME, sumV2, tree.CHLD[2].NAME, sumV3)
	if math.IsNaN(sumV1) || math.IsNaN(sumV2) || math.IsNaN(sumV3) {
		fmt.Println(tree.CHLD[0].NAME, sumV1, tree.CHLD[1].NAME, sumV2, tree.CHLD[2].NAME, sumV3)
		os.Exit(0)
	}
	tree.CHLD[0].LEN = sumV1
	tree.CHLD[1].LEN = sumV2
	tree.CHLD[2].LEN = sumV3
}
