package main

import (
	"cophymaru"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime/pprof"
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
	printFreqArg := flag.Int("pr", 10000, "Frequency with which to print to the screen")
	sampFreqArg := flag.Int("samp", 1000, "Frequency with which to sample from the chain")
	runNameArg := flag.String("o", "cophymaru", "specify the prefix for outfile names")
	threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	workersArg := flag.Int("W", 4, "Number of Go workers to use for LL calculation concurrency")
	clustArg := flag.Float64("a", 1.0, "clumpiness parameter for trait clustering algorithm")
	flag.Parse()
	f, err := os.Create("profile.prof")
	if err != nil {
		log.Fatal(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()
	//var ntax,ntraits int
	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	traits, ntax, ntraits := cophymaru.ReadContinuous(*traitArg)
	fmt.Println("SUCCESSFULLY READ IN ALIGNMENT CONTAINING ", ntax, "TAXA")
	cophymaru.MapContinuous(tree, traits, ntraits)
	rand.Seed(time.Now().UTC().UnixNano())
	for _, n := range tree.PreorderArray()[1:] {
		r := rand.Float64()
		n.LEN = r
	}
	cophymaru.IterateBMLengths(tree, *iterArg) //going to optimize branch lengths to set mean parameter for tree length in dirichlet prior
	treeOutFile := *runNameArg
	treeOutFile += ".t"
	logOutFile := *runNameArg
	logOutFile += ".mcmc"

	/*
		var mult bool
		if *threadArg == 1 {
			mult = false
		} else if *threadArg > 1 {
			mult = true
		} else {
			fmt.Println("Please pick a valid number of cores to use for the run.")
			os.Exit(0)
		}
	*/
	vcv := cophymaru.SetIdentityMatrix(ntax)
	mu0 := cophymaru.GIWStartingSampleMean(traits)
	brownianPrior := cophymaru.InitVCVPrior(mu0, 1, vcv, float64(ntax)) // instantiate the prior on the brownian proces
	chain := cophymaru.InitGibbs(traits, brownianPrior, *genArg, *printFreqArg, *sampFreqArg, *threadArg, *workersArg, *clustArg)
	//os.Exit(0)
	start := time.Now()
	chain.Run()
	elapsed := time.Since(start)
	fmt.Println("COMPLETED ", *genArg, "MCMC SIMULATIONS IN ", elapsed)
}
