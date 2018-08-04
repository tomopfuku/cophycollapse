package main

import (
	"cophycollapse"
	"flag"
	"fmt"
	"math/rand"
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
	kArg := flag.Int("K", 2, "number of clusters")
	printFreqArg := flag.Int("pr", 100, "Frequency with which to print to the screen")
	searchArg := flag.Int("f", 0, "Indicate whether to perform:\n 0\tExpectation-Maximization \nor\n1\tGreedy Hill-Climbing clustering algorithm\n2\tInitialize greedy hill climb with an EM clustering\n3\tStart greedy hill climb from randomly assigned state wth starting dimensionality K")
	//sampFreqArg := flag.Int("samp", 1, "Frequency with which to sample from the chain")
	runNameArg := flag.String("o", "cophycollapse", "specify the prefix for outfile names")
	critArg := flag.Int("c", 0, "Criterion to use for hill climbing:\n0\tAIC\n1\tBIC\n2\tAICc\n")
	//threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	//workersArg := flag.Int("W", 4, "Number of Go workers to use for LL calculation concurrency")
	//clustArg := flag.Float64("a", 1.0, "clumpiness parameter for trait clustering algorithm")
	flag.Parse()
	/*f, err := os.Create("profile.prof")
	if err != nil {
		log.Fatal(err)
	}
	//pprof.StartCPUProfile(f)
	//defer pprof.StopCPUProfile()*/
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
	cophycollapse.InitMissingValues(tree.PreorderArray())
	cophycollapse.MissingTraitsEM(tree, 100) //going to optimize branch lengths to set mean parameter for tree length in dirichlet prior
	fmt.Println(tree.Newick(true))

	fmt.Println("Starting tree AIC/BIC:", cophycollapse.CalcTreeAIC(tree, *critArg))
	//fmt.Println(tree.Newick(true))
	treeOutFile := *runNameArg
	//treeOutFile += ".trees"
	//logOutFile := *runNameArg
	//logOutFile += ".mcmc"
	if *searchArg == 0 {
		search := cophycollapse.InitEMSearch(tree, *genArg, *kArg, *printFreqArg)
		fmt.Println(search.ClusterString())
		start := time.Now()
		search.Run()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
	} else if *searchArg == 1 {
		search := cophycollapse.InitGreedyHC(tree, *genArg, *printFreqArg, *critArg, false, *kArg, treeOutFile)
		fmt.Println(search.ClusterString())
		start := time.Now()
		search.PerturbedRun()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
	} else if *searchArg == 2 {
		search := cophycollapse.InitEMSearch(tree, *genArg, *kArg, *printFreqArg)
		fmt.Println(search.ClusterString())
		start := time.Now()
		search.Run()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
		ghc := cophycollapse.TransferGreedyHC(search.Tree, search.Gen, search.PrintFreq, *critArg, search.Clusters, search.SiteAssignments, treeOutFile)
		start = time.Now()
		ghc.PerturbedRun()
		elapsed = time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)

	} else if *searchArg == 3 {
		search := cophycollapse.InitGreedyHC(tree, *genArg, *printFreqArg, *critArg, true, *kArg, treeOutFile)
		fmt.Println(search.ClusterString())
		start := time.Now()
		search.PerturbedRun()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
	}
}
