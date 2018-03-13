package main

import (
	"cophymaru"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
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
	genArg := flag.Int("gen", 500000, "number of MCMC generations to run")
	brPrior := flag.String("bl", "0", "\tspecifies the prior to place on branch lengths. \n\tOptions:\n\n\t\t0    Flat\n\t\t1    Exponential (mean = 10)\n\t\t2    Compound Dirichlet")
	fosArg := flag.String("fos", "", "file containing names of fossil tips\n\ttips should be comma-separated")
	startArg := flag.String("st", "0", "\tSpecify whether to use ML or random starting tree and branch lengths\n\t\t0    ML\n\t\t1    random")
	missingArg := flag.Bool("mis", true, "indicate whether the character matrix contains missing sites.")
	printFreqArg := flag.Int("pr", 10000, "Frequency with which to print to the screen")
	sampFreqArg := flag.Int("samp", 1000, "Frequency with which to sample from the chain")
	weightLLArg := flag.String("w", "flat", "indicate whether to place fossils using the weighted LL algorithm\nOptions:\nflat (default): flat weights\nfloat: weight site log-likelihoods using a multiplier between 0 and 1\nint: weight site log-likelihoods using a multiplier between 0 and 100\n\n")
	runNameArg := flag.String("o", "cophymaru", "specify the prefix for outfile names")
	threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	workersArg := flag.Int("W", 4, "Number of Go workers to use for LL calculation concurrency")
	algArg := flag.String("f", "0", "indicate which MCMC algorithm to perform:\n0:\tfossil placement\n1:\tsingle partition branch length estimation\n2\tcluster traits by branch length variance\n")
	clustArg := flag.Float64("a", 1.0, "clumpiness parameter for trait clustering algorithm")
	//blMeanArg := flag.Float64("beta", 1.0, "mean branch length for exponential prior or tree length for Dirichlet prior")
	flag.Parse()
	f, err := os.Create("profile")
	if err != nil {
		log.Fatal(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()
	//var ntax,ntraits int
	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	//fmt.Println("SUCCESSFULLY READ IN TREE ", tree.Newick(true))
	traits, ntax, ntraits := cophymaru.ReadContinuous(*traitArg)
	fmt.Println("SUCCESSFULLY READ IN ALIGNMENT CONTAINING ", ntax, "TAXA")
	cophymaru.MapContinuous(tree, traits, ntraits)
	//cophymaru.IterateBMLengths(tree, *iterArg)

	var weights []float64
	var fosSlice []string // read in fossil names from command line
	if *algArg == "0" {
		cophymaru.MissingTraitsEM(tree, *iterArg)
		if *weightLLArg != "flat" {
			fmt.Println("Calibrating weights to filter for concordant sites...")
			weights = cophymaru.CalibrateSiteWeights(tree, *weightLLArg)
			fmt.Println("Generating starting tree with weights:", weights)
		} else {
			for i := 0; i < ntraits; i++ {
				weights = append(weights, 1.0)
			}
		}
		for _, i := range strings.Split(*fosArg, ",") {
			fosSlice = append(fosSlice, i)
		}
		if *startArg == "0" {
			starttr, startll := cophymaru.InsertFossilTaxa(tree, traits, fosSlice, *iterArg, *missingArg, weights)
			fmt.Println("STARTING ML TREE:\n", starttr, "\n\nSTARTING MCMC WITH LOG-LIKELIHOOD ", startll)
			//fmt.Println(STARTING MCMC WITH LOG-LIKELIHOOD ", startll)
		} else if *startArg == "1" {
			cophymaru.MakeRandomStartingBranchLengths(tree)
			cophymaru.InsertFossilTaxaRandom(tree, traits, fosSlice, *iterArg, *missingArg)
			fmt.Println("\n\nINITIATING MCMC WITH RANDOM STARTING TREE") //", tree.Newick(true))
		}
	} else if *algArg == "1" {
		if *fosArg != "" {
			fmt.Println("cannot specify fossil tips when estimating branch lengths on a fixed tree!")
			os.Exit(0)
		}
		for i := 0; i < ntraits; i++ {
			weights = append(weights, 1.0) //use flat weights
		}
		cophymaru.MakeRandomStartingBranchLengths(tree)
		cophymaru.IterateBMLengths(tree, *iterArg)
	} else if *algArg == "2" {
		cophymaru.MakeRandomStartingBranchLengths(tree)

	}
	//l1 := cophymaru.CalcUnrootedLogLike(tree, true)
	//l2 := cophymaru.WeightedUnrootedLogLike(tree, true, weights)
	//fmt.Println("START COMPARISON:", l1, l2)
	//TODO: need to init new MCMC class to make this work again.
	treeOutFile := *runNameArg
	treeOutFile += ".t"
	logOutFile := *runNameArg
	logOutFile += ".mcmc"

	var mult bool
	if *threadArg == 1 {
		mult = false
	} else if *threadArg > 1 {
		mult = true
	} else {
		fmt.Println("Please pick a valid number of cores to use for the run.")
		os.Exit(0)
	}
	chain := cophymaru.InitMCMC(*genArg, treeOutFile, logOutFile, *brPrior, *printFreqArg, *sampFreqArg, *threadArg, *workersArg, mult, weights, tree, fosSlice, *algArg)
	if *algArg == "2" { // need to perform some extra steps when using the clustering algorithm
		cophymaru.InitializeClusters(chain)
		chain.NSITES = float64(ntraits)
		//chain.CLUS = c
		chain.ALPHA = *clustArg
		denom := chain.NSITES - 1 + chain.ALPHA
		chain.ALPHAPROB = (chain.ALPHA / 2.) / denom
	}
	start := time.Now()
	chain.Run()
	elapsed := time.Since(start)
	fmt.Println("COMPLETED ", *genArg, "MCMC SIMULATIONS IN ", elapsed)
}
