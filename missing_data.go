package cophycollapse

//InitMissingValues will find the missing sites in a data matrix and plug in values corresponding to the mean of the remaining sites
func InitMissingValues(tree []*Node) {
	means := CalcSiteMeans(tree)
	for _, n := range tree {
		if len(n.CHLD) == 0 {
			MakeMissingMeansTip(n, means)
		}
	}
}

//MakeMissingMeansTip will replace missing values with the mean across all tips for a single tip
func MakeMissingMeansTip(n *Node, means []float64) {
	for i := range n.CONTRT {
		if n.MIS[i] == true {
			n.CONTRT[i] = means[i]
		}
	}
}

//CalcSiteMeans will calculate the mean value for all the sites in the matrix for which the site is not missing
func CalcSiteMeans(nodes []*Node) (siteSum []float64) {
	var ntraits []int
	for range nodes[0].CONTRT {
		siteSum = append(siteSum, 0.0)
		ntraits = append(ntraits, 0)
	}
	for _, n := range nodes {
		if len(n.CHLD) != 0 {
			continue
		}
		for i, tr := range n.CONTRT {
			if n.MIS[i] != true {
				siteSum[i] += tr
				ntraits[i]++
			}
		}
	}
	for i := range siteSum {
		siteSum[i] = siteSum[i] / float64(ntraits[i])
	}
	return
}

//CalcExpectedTraitsSub will plug in the expected values for missing traits under BM using the pruning/PIC ancestral state estimation approach
func CalcExpectedTraitsSub(tree *Node, sites []int) {
	rnodes := tree.PreorderArray()
	lnode := 0
	var expect float64
	var bot float64
	var top float64
	var ltrait float64
	var llen float64
	for ind, newroot := range rnodes { //visit each internal node and check any leaves for missing data. if found, calculate input expected value as the PIC at newroot
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.CHLD {
			if len(cn.CHLD) == 0 { // visit any leaves subtending from newroot
				for traitIndex := range sites {
					if cn.MIS[traitIndex] == true { // check if trait is missing and calculate expectation for each missing value
						bot = 0.
						childCount := 0
						for _, cn2 := range tree.CHLD { //calculate PIC at newroot
							if cn2 != cn {
								childCount++
								BMPruneRootedSingle(cn2, traitIndex) // prune root to 3-tip tree
								if cn2.PRNLEN == 0. {
									cn2.PRNLEN = 0.0001
								}
								bot += 1. / cn2.PRNLEN
								if childCount == 2 {
									top = ((1. / cn2.PRNLEN) * ltrait) + ((1. / llen) * cn2.CONTRT[traitIndex])
									break
								}
								ltrait = cn2.CONTRT[traitIndex]
								llen = cn2.PRNLEN
							}
						}
						expect = top / bot
						cn.CONTRT[traitIndex] = expect
					}
				}
			}
		}
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}

//CalcExpectedTraits will plug in the expected values for missing traits under BM using the pruning/PIC ancestral state estimation approach
func CalcExpectedTraits(tree *Node) {
	rnodes := tree.PreorderArray()
	lnode := 0
	var expect float64
	var bot float64
	var top float64
	var ltrait float64
	var llen float64
	for ind, newroot := range rnodes { //visit each internal node and check any leaves for missing data. if found, calculate input expected value as the PIC at newroot
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.CHLD {
			if len(cn.CHLD) == 0 { // visit any leaves subtending from newroot
				for traitIndex := range cn.CONTRT {
					if cn.MIS[traitIndex] == true { // check if trait is missing and calculate expectation for each missing value
						bot = 0.
						childCount := 0
						for _, cn2 := range tree.CHLD { //calculate PIC at newroot
							if cn2 != cn {
								childCount++
								BMPruneRootedSingle(cn2, traitIndex) // prune root to 3-tip tree
								bot += 1. / cn2.PRNLEN
								if childCount == 2 {
									top = ((1. / cn2.PRNLEN) * ltrait) + ((1. / llen) * cn2.CONTRT[traitIndex])
									break
								}
								ltrait = cn2.CONTRT[traitIndex]
								llen = cn2.PRNLEN
							}
						}
						expect = top / bot
						cn.CONTRT[traitIndex] = expect
					}
				}
			}
		}
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}
