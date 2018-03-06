package cophymaru

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"time"
)

func brlenSlidingWindow(theta float64, wsize float64) (thetaStar float64) {
	u := rand.Float64()
	thetaStar = theta - (wsize / 2.) + (wsize * u)
	if thetaStar < 0. {
		thetaStar = -thetaStar
	}
	return
}

func cladeBrlenMultiplierProp(theta []float64, epsilon float64) (thetaStar []float64, propRat float64) {
	u := rand.Float64()
	//epsilon := 0.2
	c := math.Exp(((u - 0.5) * epsilon))
	m := float64(len(theta))
	propRat = math.Pow(c, m)
	var tmpstar float64
	for i := range theta {
		tmpstar = theta[i] * c
		thetaStar = append(thetaStar, tmpstar)
	}
	return
}

func singleBrlenMultiplierProp(theta float64, epsilon float64) (thetaStar, propRat float64) {
	min := math.Log(0.001)
	u := rand.Float64()
	//epsilon := 0.2
	c := math.Exp(((u - 0.5) * epsilon))
	thetaStar = theta * c
	propRat = c
	if math.Log(thetaStar) < min { //place a lower constraint on brlen
		thetaStar = math.Exp(2*min - math.Log(thetaStar))
		propRat = thetaStar / theta
	}
	return
}

func getProposedBrlens(nodes []*Node) []float64 {
	var oldL []float64
	for _, i := range nodes {
		oldL = append(oldL, i.LEN)
		i.LEN = brlenSlidingWindow(i.LEN, 0.2)
	}
	return oldL
}

func getFossilNodesFromLabel(fnames []string, nodes []*Node) []*Node {
	var ret []*Node
	for _, n := range nodes {
		for _, lab := range fnames {
			if n.NAME == lab {
				ret = append(ret, n)
			}
		}
	}
	return ret
}

func writeTreeFile(newick string, outFl *bufio.Writer) {
	fmt.Fprint(outFl, newick+"\n")
}

func initOutfile(flnm string) *bufio.Writer {
	f, err := os.Create(flnm)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	return w
}

//InitMCMC sets up all of the attributes of the MCMC run
func InitMCMC(gen int, treeOut, logOut string, branchPrior string, printFreq, writeFreq, maxProc, workers int, multi bool, weights []float64) (chain *MCMC) {
	chain = new(MCMC)
	chain.NGEN = gen
	chain.TREEOUTFILE = treeOut
	chain.LOGOUTFILE = logOut
	chain.BRANCHPRIOR = branchPrior
	chain.PRINTFREQ = printFreq
	chain.WRITEFREQ = writeFreq
	chain.PROC = maxProc
	chain.WORKERS = workers
	chain.MULTI = multi
	chain.SITEWEIGHTS = weights
	return
}

//initializeRun sets global variables for the run
func (chain *MCMC) initializeRun(tree *Node, fnames []string) (nodes, fos, innodes []*Node, prior *BranchLengthPrior, treeLogLikelihood *LL) {
	nodes = tree.PreorderArray()[1:]
	fos = getFossilNodesFromLabel(fnames, nodes)
	innodes = InternalNodeSlice(nodes)
	InitParallelPRNLEN(nodes)
	prior = InitializePrior(chain.BRANCHPRIOR, nodes)
	treeLogLikelihood = InitLL(chain.MULTI, chain.WORKERS, chain.SITEWEIGHTS)
	runtime.GOMAXPROCS(chain.PROC)
	return
}

//MCMC is a struct for storing information about the current run
type MCMC struct {
	NGEN        int
	TREEOUTFILE string
	LOGOUTFILE  string
	BRANCHPRIOR string
	PRINTFREQ   int
	WRITEFREQ   int
	PROC        int
	WORKERS     int
	MULTI       bool
	SITEWEIGHTS []float64
}

//Run will run Markov Chain Monte Carlo simulations, adjusting branch lengths and fossil placements
func (chain *MCMC) Run(tree *Node, fnames []string) {
	f, err := os.Create(chain.TREEOUTFILE)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)
	logFile, err := os.Create(chain.LOGOUTFILE)
	if err != nil {
		log.Fatal(err)
	}
	logWriter := bufio.NewWriter(logFile)
	nodes, fos, inNodes, branchPrior, treeLogLikelihood := chain.initializeRun(tree, fnames)
	ll := treeLogLikelihood.Calc(tree, true)
	lp := branchPrior.Calc(nodes)
	acceptanceCount := 0.0
	topAcceptanceCount := 0.
	var oldLL, acceptanceRatio, topAcceptanceRatio float64
	epsilon := 0.1
	for i := 0; i < chain.NGEN; i++ {
		oldLL = ll
		if i%3 == 0 || i == 0 {
			//lp, ll = singleBranchLengthUpdate(ll, lp, nodes, inNodes, tree, branchPrior, missing, weights) //NOTE uncomment to sample BRLENS
			lp, ll = fossilPlacementUpdate(ll, lp, fos, nodes, tree, chain, branchPrior, treeLogLikelihood)
			if ll != oldLL {
				topAcceptanceCount = topAcceptanceCount + 1.0
			}
		} else {
			s1 := rand.NewSource(time.Now().UnixNano())
			r1 := rand.New(s1)
			r := r1.Float64()
			if r > 0.1 { // apply single branch length update 95% of the time
				lp, ll = singleBranchLengthUpdate(ll, lp, nodes, inNodes, tree, chain, branchPrior, treeLogLikelihood, epsilon) //NOTE uncomment to sample BRLENS
			} else {
				lp, ll = cladeBranchLengthUpdate(ll, lp, nodes, inNodes, tree, chain, branchPrior, treeLogLikelihood, epsilon)
			}
			if ll != oldLL {
				acceptanceCount = acceptanceCount + 1.0
			}
		}

		if i%200 == 0 && i <= 10000 && i != 0 { // use burn in period to adjust the branch length multiplier step length every 200 generations
			acceptanceRatio = acceptanceCount / float64(i)
			epsilon = adjustBranchLengthStepLength(epsilon, acceptanceRatio)
		}
		if i == 0 {
			fmt.Println("generation", "logPrior", "logLikelihood", "acceptanceRatio", "topologyAcceptanceRatio")
			fmt.Println("0", lp, ll, "NA", "NA")
		}
		if i%chain.PRINTFREQ == 0 && i != 0 {
			acceptanceRatio = acceptanceCount / float64(i)
			topAcceptanceRatio = topAcceptanceCount / float64(i)
			fmt.Println(i, lp, ll, acceptanceRatio, topAcceptanceRatio)
		}

		if i%chain.WRITEFREQ == 0 {
			fmt.Fprint(logWriter, strconv.Itoa(i)+"\t"+strconv.FormatFloat(lp, 'f', -1, 64)+"\t"+strconv.FormatFloat(ll, 'f', -1, 64)+"\n")
			/*
				for _, ln := range nodes {
					fmt.Fprint(lw, strconv.FormatFloat(ln.LEN, 'f', -1, 64)+"\t")
				}
				fmt.Fprint(lw, "\n")
			*/
			//writeTreeFile(tree.Newick(true),w)
			fmt.Fprint(w, tree.Newick(true)+";\n")
		}
	}
	logWriter.Flush()
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

func adjustBranchLengthStepLength(epsilon, acceptanceRatio float64) (epsilonStar float64) { //this will calculate the optimal step length for the single branch length multiplier proposal
	acceptanceRatioStar := 0.44 // this is the optimal acceptance probability for uniform proposals
	s := math.Pi / 2.
	epsilonStar = epsilon * (math.Tan(s*acceptanceRatio) / math.Tan(s*acceptanceRatioStar))
	return
}

//this move prunes and regrafts a fossil, creating random variables for the altered branch lengths
//the procedure is basically the same as the SPR move as described in the Yang (2014) Mol. Evol. book
func fossilPlacementUpdate(ll, lp float64, fnodes, nodes []*Node, tree *Node, chain *MCMC, branchPrior *BranchLengthPrior, treeLL *LL) (float64, float64) {
	fn := drawRandomNode(fnodes) //draw a random fossil tip
	reattach := drawRandomReattachment(fn, nodes)
	x, p, lastn := PruneFossilTip(fn)
	r := GraftFossilTip(fn.PAR, reattach)
	propRat := r / (x * p)
	//tree.UnmarkAll()
	//lastn.UnmarkToRoot(tree)
	//fn.UnmarkToRoot(tree)
	//reattach.UnmarkToRoot(tree)
	//llstar := WeightedUnrootedLogLike(tree, true, weights)
	llstar := treeLL.Calc(tree, true)
	//llstar1 := WeightedUnrootedLogLike(tree, true, weights)
	//MarkAll(nodes)
	//fmt.Println(llstar, llstar1)
	lpstar := branchPrior.Calc(nodes)
	alpha := math.Exp(lpstar-lp) * math.Exp(llstar-ll) * propRat
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r = r1.Float64()
	if r < alpha {
		ll = llstar
		lp = lpstar
	} else { //move fossil back to its previous position and restore old branch lengths
		PruneFossilTip(fn)
		GraftFossilTip(fn.PAR, lastn)
		lastn.LEN = x
		fn.PAR.LEN = p
		//fn.UnmarkToRoot(tree)
		//reattach.UnmarkToRoot(tree)
		//tree.UnmarkAll()
	}
	return lp, ll
}

// here is a test
func singleBranchLengthUpdate(ll, lp float64, nodes, inNodes []*Node, tree *Node, chain *MCMC, branchPrior *BranchLengthPrior, treeLL *LL, epsilon float64) (float64, float64) {
	updateNode := RandomNode(nodes)
	soldL := updateNode.LEN
	//updateNode.UnmarkToRoot(tree)
	var propRat float64
	updateNode.LEN, propRat = singleBrlenMultiplierProp(updateNode.LEN, epsilon)
	//updateNode.UnmarkToRoot(tree)
	llstar := treeLL.Calc(tree, true)
	//MarkAll(nodes)

	lpstar := branchPrior.Calc(nodes)

	alpha := math.Exp(lpstar-lp) * math.Exp(llstar-ll) * propRat
	//fmt.Println(llstar, ll, llstar-ll)
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r := r1.Float64()
	if r < alpha {
		ll = llstar
		lp = lpstar
	} else {
		//updateNode.UnmarkToRoot(tree)
		updateNode.LEN = soldL
	}
	return lp, ll
}

func cladeBranchLengthUpdate(ll, lp float64, nodes, inNodes []*Node, tree *Node, chain *MCMC, branchPrior *BranchLengthPrior, treeLL *LL, epsilon float64) (float64, float64) {
	updateNode := RandomInternalNode(nodes)
	updateClade := updateNode.PostorderArray()
	var oldlens []float64
	for _, node := range updateClade {
		oldlens = append(oldlens, node.LEN)
	}
	newlens, propRat := cladeBrlenMultiplierProp(oldlens, epsilon)
	for i, node := range updateClade {
		node.LEN = newlens[i]
		node.MRK = false
	}
	//updateNode.UnmarkToRoot(tree)
	llstar := treeLL.Calc(tree, true)
	lpstar := branchPrior.Calc(nodes)

	alpha := math.Exp(lpstar-lp) * math.Exp(llstar-ll) * propRat
	//fmt.Println(llstar, ll, llstar-ll)
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r := r1.Float64()
	if r < alpha {
		ll = llstar
		lp = lpstar
	} else {
		//updateNode.UnmarkToRoot(tree)
		for i, node := range updateClade {
			node.LEN = oldlens[i]
		}
	}
	return lp, ll
}

func drawRandomNode(n []*Node) (rnode *Node) {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rnoden := r1.Intn(len(n))
	rnode = n[rnoden]
	return
}

//RandomNode will pull a random node from a slice of nodes
func RandomNode(nodes []*Node) (rnode *Node) {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rnoden := r1.Intn(len(nodes))
	rnode = nodes[rnoden] //choose a random node
	return
}

//RandomInternalNode will draw a random internal node
func RandomInternalNode(nodes []*Node) (rnode *Node) {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rnoden := r1.Intn(len(nodes))
	rnode = nodes[rnoden] //choose a random reattachment point
	if len(rnode.CHLD) == 0 {
		rnode = RandomInternalNode(nodes)
	}
	return
}

func drawRandomReattachment(fn *Node, nodes []*Node) (rnode *Node) {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rnoden := r1.Intn(len(nodes))
	rnode = nodes[rnoden] //choose a random reattachment point
	if rnode == fn {
		rnode = drawRandomReattachment(fn, nodes)
	} else if rnode == fn.PAR {
		rnode = drawRandomReattachment(fn, nodes)
	}
	return
}

//TODO: should probably combine these two into a single SPR move function

//PruneFossilTip removes a fossil tip, along with its parent node from a tree, returning the new branch length left by the gap
//return value is for calculating proposal ratio later
func PruneFossilTip(tip *Node) (x, p float64, lastn *Node) {
	var n *Node
	newpar := tip.PAR
	for _, chd := range newpar.CHLD {
		if chd != tip {
			n = chd
		}
	}
	newpar.PAR.AddChild(n)
	newpar.RemoveChild(n)
	n.PAR = newpar.PAR
	newpar.PAR.RemoveChild(newpar)
	x = n.LEN
	p = newpar.LEN
	lastn = n
	n.LEN = n.LEN + newpar.LEN
	return
}

//GraftFossilTip reattaches a fossil tip (newpar) to branch n at a random point
func GraftFossilTip(newpar *Node, n *Node) float64 {
	r := n.LEN
	n.PAR.AddChild(newpar)
	n.PAR.RemoveChild(n)
	newpar.AddChild(n)
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	u := r1.Float64()
	newpar.LEN = u * n.LEN
	n.LEN = n.LEN * (1 - u)
	return r
}
