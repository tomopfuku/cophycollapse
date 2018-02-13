package main

import (
	"cophymaru"
	"flag"
	"fmt"
	"strings"
	"time"
)

func postorder(curnode *cophymaru.Node) {
	for _, chld := range curnode.CHLD {
		postorder(chld)
	}
	fmt.Println(curnode.NAME, curnode.CONTRT, curnode.LEN)
}

func main() {
	treeArg := flag.String("t", "", "input tree")
	traitArg := flag.String("m", "", "continuous traits")
	iterArg := flag.Int("i", 5, "num iterations for branch length iteration")
	genArg := flag.Int("gen", 5000, "number of MCMC generations to run")
	brPrior := flag.String("blpr", "0", "\tspecifies the prior to place on branch lengths. \n\tOptions:\n\n\t\t0    Flat\n\t\t1    Exponential (mean = 10)\n")
	fosArg := flag.String("fos", "", "file containing names of fossil tips\n\ttips should be comma-separated")
	startArg := flag.String("st", "0", "\tSpecify whether to use ML or random starting tree and branch lengths\n\t\t0    ML\n\t\t1    random\n")
	missingArg := flag.Bool("mis", true, "indicate whether the character matrix contains missing sites.")
	printFreqArg := flag.Int("pr", 100, "Frequency with which to print to the screen")
	sampFreqArg := flag.Int("samp", 100, "Frequency with which to sample from the chain")
	weightLLArg := flag.Bool("w", false, "indicate whether to place fossils using the weighted LL algorithm")
	flag.Parse()
	//var ntax,ntraits int
	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	fmt.Println("READ IN TREE: ", tree.Newick(true))
	traits, ntax, ntraits := cophymaru.ReadContinuous(*traitArg)
	fmt.Println("CONTAINING ", ntax, "TAXA")
	cophymaru.MapContinuous(tree, traits, ntraits)
	/*/ test random vs reference tree LL comparison
	randTree := cophymaru.RandomUnrootedTree(tree)
	cophymaru.IterateBMLengths(tree, *iterArg)
	cophymaru.IterateBMLengths(randTree, *iterArg)
	l1 := cophymaru.MissingUnrootedLogLike(randTree, true)
	l2 := cophymaru.MissingUnrootedLogLike(tree, true)
	fmt.Println(randTree.Newick(true))
	fmt.Println(l1, l2)
	os.Exit(0)
	*/
	cophymaru.IterateBMLengths(tree, *iterArg)
	var weights []float64
	if *weightLLArg == true {
		fmt.Println("Calibrating weights to filter for concordant sites...")
		weights = cophymaru.CalibrateSiteWeights(tree)
		//fmt.Println("Generating starting tree with weights:", weights)
	} else {
		for i := 0; i < ntraits; i++ {
			weights = append(weights, 1.0)
		}
	}
	var fosSlice []string // read in fossil names from command line
	for _, i := range strings.Split(*fosArg, ",") {
		fosSlice = append(fosSlice, i)
	}
	if *startArg == "0" {
		starttr, startll := cophymaru.InsertFossilTaxa(tree, traits, fosSlice, *iterArg, *missingArg, weights)
		fmt.Println("STARTING ML TREE:\n", starttr, "\n\nSTARTING MCMC WITH LOG-LIKELIHOOD ", startll)
	} else if *startArg == "1" {
		cophymaru.MakeRandomStartingBranchLengths(tree)
		cophymaru.InsertFossilTaxaRandom(tree, traits, fosSlice, *iterArg, *missingArg)
		fmt.Println("START: ", tree.Newick(true))
	}
	//l1 := cophymaru.CalcUnrootedLogLike(tree, true)
	//l2 := cophymaru.WeightedUnrootedLogLike(tree, true, weights)
	//fmt.Println("START COMPARISON:", l1, l2)
	start := time.Now()
	cophymaru.MCMC(tree, *genArg, fosSlice, "tmp/test.t", "tmp/test.mcmc", *brPrior, *missingArg, *printFreqArg, *sampFreqArg, weights)
	elapsed := time.Since(start)
	fmt.Println("COMPLETED ", *genArg, "MCMC SIMULATIONS IN ", elapsed)
}