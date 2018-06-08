package cophycollapse

//TODO: this has a bug somewhere. need to track it down at some point.

//InsertFossilTaxa finds the ML placement of a set of fossils via stepwise addition
func InsertFossilTaxa(tree *Node, traits map[string][]float64, fosNms []string, iter int, missing bool, weights []float64) (besttr string, bestll float64) {
	nodes := tree.PreorderArray()
	//TODO: get this missing data stuff up and running
	//trMeans := CalcSiteMeans(nodes)
	//MakeMissingMeansArray(nodes)
	for _, curFos := range fosNms {
		ftip := new(Node)
		ftip.NAME = curFos
		ftip.LEN = 0.01
		ftip.CONTRT = traits[curFos]
		MakeMissingDataSlice(ftip)
		//MakeMissingMeansTip(ftip, trMeans) //replace missing traits with the mean value across all sample traits
		newpar := new(Node)
		newpar.LEN = 0.01
		newpar.AddChild(ftip)
		for range ftip.CONTRT {
			newpar.CONTRT = append(newpar.CONTRT, float64(0.0))
			newpar.MIS = append(newpar.MIS, false)
			newpar.LL = append(newpar.LL, 0.)
		}
		bestll = -1000000000000.0
		besttr = ""
		curll := float64(0.0)
		var bestPlace *Node
		for _, n := range nodes[1:] {
			GraftFossilTip(newpar, n)
			//fmt.Println(tree.Newick(true))
			//start := time.Now()
			if missing == false {
				IterateBMLengths(tree, iter)
				curll = CalcUnrootedLogLike(tree, true)
			} else if missing == true {
				MissingTraitsEM(tree, iter)
				curll = WeightedUnrootedLogLike(tree, true, weights)
			}
			//end := time.Now()
			//fmt.Println(end.Sub(start))

			if curll > bestll {
				bestll = curll
				besttr = tree.Newick(true)
				bestPlace = n
			}
			//fmt.Println(curll)
			//fmt.Println(tree.Newick(true))
			PruneTip(newpar, n)
		}
		//fmt.Println(bestll,besttr)
		GraftFossilTip(newpar, bestPlace)
		if missing == false {
			IterateBMLengths(tree, iter)
		} else if missing == true {
			MissingTraitsEM(tree, iter)
		}
		//nodes = tree.PreorderArray() // need to reinitialize the node list to include the now-placed fossil
		//fmt.Println(tree.Newick(true))
	}
	return
}

//InsertFossilTaxaRandom will randomly insert all of the fossil taxa in a dataset
func InsertFossilTaxaRandom(tree *Node, traits map[string][]float64, fosNms []string, iter int, missing bool) {
	nodes := tree.PreorderArray()
	for _, curFos := range fosNms {
		ftip := new(Node)
		ftip.NAME = curFos
		ftip.LEN = 0.01
		ftip.CONTRT = traits[curFos]
		MakeMissingDataSlice(ftip)
		//MakeMissingMeansTip(ftip, trMeans) //replace missing traits with the mean value across all sample traits
		newpar := new(Node)
		newpar.LEN = 0.01
		newpar.AddChild(ftip)
		for range ftip.CONTRT {
			newpar.CONTRT = append(newpar.CONTRT, float64(0.0))
			newpar.MIS = append(newpar.MIS, false)
			newpar.LL = append(newpar.LL, 0.)
		}
		reattach := RandomNode(nodes[1:])
		GraftFossilTip(newpar, reattach)
		nodes = tree.PreorderArray() // need to reinitialize the node list to include the now-placed fossil
	}
	if missing == false {
		IterateBMLengths(tree, iter)
	} else if missing == true {
		MissingTraitsEM(tree, iter)
	}
}

//PruneTip will prune a fossil (or any) tip from a tree
//newpar is the parent of the tip you wish to prune, n is the other child of newpar
func PruneTip(newpar *Node, n *Node) {
	newpar.PAR.AddChild(n)
	newpar.RemoveChild(n)
	n.PAR = newpar.PAR
	newpar.PAR.RemoveChild(newpar)
	n.LEN = n.LEN + newpar.LEN
}

//GraftTip is the inverse of above
func GraftTip(newpar *Node, n *Node) {
	n.PAR.AddChild(newpar)
	n.PAR.RemoveChild(n)
	newpar.AddChild(n)
	newpar.LEN = n.LEN / 2.
	n.LEN = newpar.LEN
}
