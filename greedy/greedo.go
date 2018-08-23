package main

import (
	"bufio"
	"cophycollapse"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
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
	mclArg := flag.String("start", "", "user-specified starting clusters")
	genArg := flag.Int("gen", 500000, "number of MCMC generations to run")
	kArg := flag.Int("K", 2, "maximum number of clusters")
	minKArg := flag.Int("minK", 1, "minimum number of clusters")
	printFreqArg := flag.Int("pr", 100, "Frequency with which to print to the screen")
	searchArg := flag.Int("f", 3, "0\tOptimize branch lengths for a user-specified clustering\n1\tOutput distance matrices calculated for each cluster provided by the -start argument\n3\tPerform cluster analysis")
	//sampFreqArg := flag.Int("samp", 1, "Frequency with which to sample from the chain")
	runNameArg := flag.String("o", "cophycollapse", "specify the prefix for outfile names")
	critArg := flag.Int("c", 0, "Criterion to use for hill climbing:\n0\tAIC\n1\tBIC\n2\tAICc\n")
	splitGenArg := flag.Int("split", 10, "Number of iterations to run at each splitting step")
	//threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	//workersArg := flag.Int("W", 4, "Number of Go workers to use for LL calculation concurrency")
	clustArg := flag.Float64("a", 1.0, "concentration parameter for new cluster penalty")
	flag.Parse()
	f, err := os.Create("profile.prof")
	if err != nil {
		log.Fatal(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()
	//var ntax,ntraits int
	runtime.GOMAXPROCS(2)
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
	//fmt.Println(tree.Newick(true))
	cophycollapse.InitParallelPRNLEN(tree.PreorderArray())
	//fmt.Println("Starting tree AIC/BIC:", cophycollapse.CalcTreeAIC(tree, *critArg))
	//fmt.Println(tree.Newick(true))
	treeOutFile := *runNameArg
	if *searchArg == 3 {
		search := cophycollapse.InitGreedyHC(tree, *genArg, *printFreqArg, *critArg, true, *kArg, treeOutFile, *splitGenArg, *clustArg, *minKArg)
		//fmt.Println(search.ClusterString())
		start := time.Now()
		search.PerturbedRun()
		elapsed := time.Since(start)
		fmt.Println("COMPLETED ", *genArg, "ITERATIONS IN ", elapsed)
	} else if *searchArg == 0 {
		if *mclArg == "" {
			fmt.Println("You need to specify a cluster input file to run this option")
			os.Exit(1)
		}
		f, err := os.Create(*runNameArg + "_TREES")
		if err != nil {
			log.Fatal(err)
		}
		w := bufio.NewWriter(f)
		clusters := cophycollapse.ReadMCLoutput(*mclArg)
		nodes := tree.PreorderArray()
		for _, c := range clusters {
			for _, n := range nodes[1:] {
				r := rand.Float64()
				n.LEN = r
			}
			cophycollapse.ClusterMissingTraitsEM(tree, c, 100)
			sites := ""
			for _, site := range c.Sites {
				sites += strconv.Itoa(site) + "\t"
			}
			fmt.Fprint(w, sites+"\n"+tree.Newick(true)+"\n")
		}
		err = w.Flush()
		if err != nil {
			log.Fatal(err)
		}
		f.Close()
	} else if *searchArg == 1 {
		if *mclArg == "" {
			fmt.Println("You need to specify a cluster input file to run this option")
			os.Exit(1)
		}

		clusters := cophycollapse.ReadMCLoutput(*mclArg)
		nodes := tree.PreorderArray()
		clmap, err := os.Create("cluster_key")
		if err != nil {
			log.Fatal(err)
		}
		w1 := bufio.NewWriter(clmap)
		for i, c := range clusters {
			f, err := os.Create(strconv.Itoa(i) + ".phy")
			if err != nil {
				log.Fatal(err)
			}
			w := bufio.NewWriter(f)
			dm := cophycollapse.SubDM(nodes, c)
			out := cophycollapse.DMtoPhylip(dm, nodes)
			fmt.Fprint(w, strings.Join(out, "\n"))
			err = w.Flush()
			if err != nil {
				log.Fatal(err)
			}
			f.Close()
			var strsites []string
			for _, i := range c.Sites {
				strsites = append(strsites, strconv.Itoa(i))
			}
			fmt.Fprint(w1, "CLUSTER"+strconv.Itoa(i)+"\t"+strings.Join(strsites, "\t")+"\n")
		}
		err = w1.Flush()
		if err != nil {
			log.Fatal(err)
		}
		clmap.Close()
	}
}
