package main

import (
	"cophycollapse"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime/pprof"
	"time"
)

func postorder(curnode *cophycollapse.Node) {
	for _, chld := range curnode.CHLD {
		postorder(chld)
	}
	fmt.Println(curnode.NAME, curnode.CONTRT, curnode.LEN)
}

func main() {
	treeArg := flag.String("t", "", "input tree")
	traitArg := flag.String("m", "", "continuous traits")
	genArg := flag.Int("gen", 500000, "number of MCMC generations to run")
	printFreqArg := flag.Int("pr", 50, "Frequency with which to print to the screen")
	sampFreqArg := flag.Int("samp", 1, "Frequency with which to sample from the chain")
	runNameArg := flag.String("o", "cophycollapse", "specify the prefix for outfile names")
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
	nwk := cophycollapse.ReadLine(*treeArg)[0]
	tree := cophycollapse.ReadTree(nwk)
	traits, ntax, ntraits := cophycollapse.ReadContinuous(*traitArg)
	fmt.Println("SUCCESSFULLY READ IN ALIGNMENT CONTAINING ", ntax, "TAXA")
	cophycollapse.MapContinuous(tree, traits, ntraits)
	rand.Seed(time.Now().UTC().UnixNano())
	for _, n := range tree.PreorderArray()[1:] {
		r := rand.Float64()
		n.LEN = r
	}
	//	cophycollapse.IterateBMLengths(tree, *iterArg) //going to optimize branch lengths to set mean parameter for tree length in dirichlet prior
	treeOutFile := *runNameArg
	treeOutFile += ".clusters"
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
	//vcv := cophycollapse.SetIdentityMatrix(ntax)
	nodes := tree.PreorderArray()
	dist := cophycollapse.DM(nodes)
	fmt.Println(dist.MatSites[0])
	ngprior := cophycollapse.InitNGPrior(0.0, 2., 1.25, 2.)
	chain := cophycollapse.InitUVNGibbs(nodes, ngprior, *genArg, *printFreqArg, *sampFreqArg, *threadArg, *workersArg, *clustArg, treeOutFile, logOutFile)
	start := time.Now()
	chain.Run()
	elapsed := time.Since(start)
	fmt.Println("COMPLETED ", *genArg, "COLLAPSED GIBBS SIMULATIONS IN ", elapsed)
}
