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
func InitMCMC(gen int, treeOut, logOut string, branchPrior string, printFreq, writeFreq, maxProc, workers int, multi bool, weights []float64, tree *Node, fnames []string, alg string) (chain *MCMC) {
	chain = new(MCMC)
	chain.TREE = tree
	chain.NODES = tree.PreorderArray()
	InitParallelPRNLEN(chain.NODES)
	chain.INNODES = InternalNodeSlice(chain.NODES)
	chain.FOS = getFossilNodesFromLabel(fnames, chain.NODES)
	chain.NGEN = gen
	chain.TREEOUTFILE = treeOut
	chain.LOGOUTFILE = logOut
	chain.BRANCHPRIOR = InitializePrior(branchPrior, chain.NODES)
	chain.PRINTFREQ = printFreq
	chain.WRITEFREQ = writeFreq
	chain.PROC = maxProc
	chain.WORKERS = workers
	chain.MULTI = multi
	chain.SITEWEIGHTS = weights
	chain.TREELL = InitLL(chain.MULTI, chain.WORKERS, chain.SITEWEIGHTS)
	chain.STEPLEN = 0.1
	chain.ALG = alg
	return
}

//MCMC is a struct for storing information about the current run
type MCMC struct {
	NGEN        int
	TREEOUTFILE string
	LOGOUTFILE  string
	BRANCHPRIOR *BranchLengthPrior
	TREELL      *LL
	PRINTFREQ   int
	WRITEFREQ   int
	PROC        int
	WORKERS     int
	MULTI       bool
	SITEWEIGHTS []float64
	NODES       []*Node
	FOS         []*Node
	INNODES     []*Node
	TREE        *Node
	STEPLEN     float64
	ALG         string
	CLUS        []int // stores the current allocation of clusters for each site.
	NSITES      float64
	ALPHA       float64
	ALPHAPROB   float64
	CLUSTERSET  map[int][]int // stores the single set of clusters and their assignments
	UNIQUEK     []int
}

//Run will run Markov Chain Monte Carlo simulations, adjusting branch lengths and fossil placements
func (chain *MCMC) Run() {
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
	runtime.GOMAXPROCS(chain.PROC)
	chain.TREELL.CUR = chain.TREELL.Calc(chain.TREE, true)
	chain.BRANCHPRIOR.CUR = chain.BRANCHPRIOR.Calc(chain.NODES)
	acceptanceCount := 0.0
	topAcceptanceCount := 0.
	var acceptanceRatio, topAcceptanceRatio float64
	for i := 0; i < chain.NGEN; i++ {
		chain.update(i, &topAcceptanceCount, &acceptanceCount)

		if i%200 == 0 && i <= 10000 && i != 0 { // use burn in period to adjust the branch length multiplier step length every 200 generations
			acceptanceRatio = acceptanceCount / float64(i)
			chain.STEPLEN = adjustBranchLengthStepLength(chain.STEPLEN, acceptanceRatio)
		}
		if i == 0 {
			fmt.Println("generation", "logPrior", "logLikelihood", "acceptanceRatio", "topologyAcceptanceRatio")
			fmt.Println("0", chain.BRANCHPRIOR.CUR, chain.TREELL.CUR, "NA", "NA")
		}
		if i%chain.PRINTFREQ == 0 && i != 0 {
			acceptanceRatio = acceptanceCount / float64(i)
			topAcceptanceRatio = topAcceptanceCount / float64(i)
			fmt.Println(i, chain.BRANCHPRIOR.CUR, chain.TREELL.CUR, acceptanceRatio, topAcceptanceRatio, chain.STEPLEN)
		}

		if i%chain.WRITEFREQ == 0 {
			fmt.Fprint(logWriter, strconv.Itoa(i)+"\t"+strconv.FormatFloat(chain.BRANCHPRIOR.CUR, 'f', -1, 64)+"\t"+strconv.FormatFloat(chain.TREELL.CUR, 'f', -1, 64)+"\n")
			/*
				for _, ln := range nodes {
					fmt.Fprint(lw, strconv.FormatFloat(ln.LEN, 'f', -1, 64)+"\t")
				}
				fmt.Fprint(lw, "\n")
			*/
			//writeTreeFile(tree.Newick(true),w)
			fmt.Fprint(w, chain.TREE.Newick(true)+";\n")
		}
	}
	logWriter.Flush()
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

func (chain *MCMC) update(i int, topAcceptanceCount *float64, acceptanceCount *float64) {
	chain.TREELL.LAST = chain.TREELL.CUR
	if chain.ALG == "0" {
		if i%3 == 0 || i == 0 {
			//lp, ll = singleBranchLengthUpdate(ll, lp, nodes, inNodes, tree, branchPrior, missing, weights) //NOTE uncomment to sample BRLENS
			chain.fossilPlacementUpdate()
			//lp, ll = fossilPlacementUpdate(ll, lp, fos, nodes, tree, chain, branchPrior, treeLogLikelihood)
			if chain.TREELL.CUR != chain.TREELL.LAST {
				*topAcceptanceCount += 1.0
			}
		} else {
			s1 := rand.NewSource(time.Now().UnixNano())
			r1 := rand.New(s1)
			r := r1.Float64()
			if r > 0.1 { // apply single branch length update 95% of the time
				chain.singleBranchLengthUpdate()
			} else {
				chain.cladeBranchLengthUpdate()
			}
			if chain.TREELL.CUR != chain.TREELL.LAST {
				*acceptanceCount += 1.0
			}
		}
	} else if chain.ALG == "1" {
		s1 := rand.NewSource(time.Now().UnixNano())
		r1 := rand.New(s1)
		r := r1.Float64()
		if r < 0.99 { // apply single branch length update 95% of the time
			chain.singleBranchLengthUpdate()
		} else {
			chain.cladeBranchLengthUpdate()
		}
		if chain.TREELL.CUR != chain.TREELL.LAST {
			*acceptanceCount += 1.0
		}
	} else if chain.ALG == "2" {
		s1 := rand.NewSource(time.Now().UnixNano())
		r1 := rand.New(s1)
		r := r1.Float64()
		if r < 0.98 { // apply single branch length update 95% of the time
			chain.singleBranchLengthUpdateCluster()
		} else {
			chain.gibbsClusterUpdate()
		}
		if chain.TREELL.CUR != chain.TREELL.LAST {
			*acceptanceCount += 1.0
		}
	}
}
func adjustBranchLengthStepLength(epsilon, acceptanceRatio float64) (epsilonStar float64) { //this will calculate the optimal step length for the single branch length multiplier proposal
	acceptanceRatioStar := 0.44 // this is the optimal acceptance probability for uniform proposals
	s := math.Pi / 2.
	epsilonStar = epsilon * (math.Tan(s*acceptanceRatio) / math.Tan(s*acceptanceRatioStar))
	return
}

func (chain *MCMC) gibbsClusterUpdate() {
	for i := range chain.CLUS {
		chain.siteClusterUpdate(i)
	}
}

//this will update the cluster assignment of a single specified site
func (chain *MCMC) siteClusterUpdate(curSite int) {
	catMinusI := chain.CLUSTERSET //make(map[int][]int)
	alone := false
	curSiteCluster := chain.CLUS[curSite]
	if len(catMinusI[curSiteCluster]) == 1 {
		delete(catMinusI, curSiteCluster)
		alone = true
	} else {
		var newsites []int
		for _, c := range catMinusI[curSiteCluster] {
			if c != curSite {
				newsites = append(newsites, c)
			}
		}
		catMinusI[curSiteCluster] = newsites
	}
	var clusterProbs map[int]float64
	if alone == false { // if curSite belongs to a cluster shared with other data
		aux := -2
		chain.drawAuxBL(aux)
		clusterProbs = chain.clusterAssignmentProbs(catMinusI, curSite, aux)
	} else { //treat curSiteCluster as an auxilliary class if curSite occupies its own cluster
		aux := -1
		chain.drawAuxBL(aux)
		clusterProbs = chain.clusterAssignmentProbs(catMinusI, curSite, aux)
	}
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r := r1.Float64()
	cumprob := 0.
	var newcluster int
	for k := range clusterProbs {
		cumprob += clusterProbs[k]
		if cumprob > r {
			newcluster = k
		}
	}
	if newcluster < 0 { // create a new cluster K+1 if curSite was assigned to one of the aux classes
		newlens := chain.TREE.ClustLEN[newcluster]
		newcluster = Max(catMinusI) + 1
		chain.TREE.ClustLEN[newcluster] = newlens
	}
	//chain.CLUS[curSite] = newcluster
	if newcluster != curSiteCluster {
		chain.updateAssignmentVector(curSite, newcluster)
		if alone == true {
			for _, n := range chain.NODES {
				delete(n.ClustLEN, curSiteCluster) // delete any clusters that are empty
			}
			delete(chain.CLUSTERSET, curSiteCluster)
			chain.updateUniqueK(curSiteCluster)
		} else {
			chain.CLUSTERSET[curSiteCluster] = catMinusI[curSiteCluster]
		}
	}
}

func (chain *MCMC) updateUniqueK(del int) {
	var new []int
	for _, c := range chain.UNIQUEK {
		if c != del {
			new = append(new, c)
		}
	}
	chain.UNIQUEK = new
}

func (chain *MCMC) updateAssignmentVector(curSite, newcluster int) {
	var newAssignmentVector []int
	for i, c := range chain.CLUS {
		if i != curSite {
			newAssignmentVector = append(newAssignmentVector, c)
		} else {
			newAssignmentVector = append(newAssignmentVector, newcluster)
		}
	}
	chain.CLUS = newAssignmentVector
}

//this could be sped up if i kept track of both the current and last cluster assignments for each site and reused the likelihood calc for sites that haven't changed clusters
func (chain *MCMC) clusterAssignmentProbs(cat map[int][]int, cur, aux int) (prob map[int]float64) {
	prob = make(map[int]float64)
	var rat float64
	ratsum := 0.
	denom := chain.NSITES - 1 + chain.ALPHA
	for k := range cat { // calculate assignment probs for all assigned categories
		rat = float64(len(cat[k])) / denom
		siteLL := SingleSiteLikeCluster(chain, cur, k)
		rat = rat * siteLL
		prob[k] = rat
		ratsum += rat
	}
	//now need to calculate the probabilites for reassigning to an auxilliary category
	rat = chain.ALPHAPROB
	var siteLL float64
	curSiteCluster := chain.CLUS[cur]
	if aux == -1 { // treat curSiteCluster as one of the auxiliary clusters if curSite has its own cluster (is alone)
		siteLL = SingleSiteLikeCluster(chain, cur, curSiteCluster)
		rat = rat * siteLL
		prob[curSiteCluster] = rat
		ratsum += rat
	}
	for k := -1; k <= aux; k-- {
		siteLL = SingleSiteLikeCluster(chain, cur, k)
		rat = rat * siteLL
		prob[k] = rat
		ratsum += rat
	}
	for k := range prob {
		prob[k] = prob[k] / ratsum
	}
	return
}

func (chain *MCMC) drawAuxBL(aux int) {
	for _, n := range chain.NODES {
		for i := -1; i <= aux; i-- {
			s1 := rand.NewSource(time.Now().UnixNano())
			r1 := rand.New(s1)
			u := r1.Float64()
			n.ClustLEN[i] = u //ClustLEN is a map, not list
		}
	}
}

func (chain *MCMC) singleBranchLengthUpdateCluster() {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	cluster := chain.UNIQUEK[r1.Intn(len(chain.UNIQUEK))]
	for _, n := range chain.NODES {
		n.LEN = n.ClustLEN[cluster]
	}
	updateNode := RandomNode(chain.NODES)
	soldL := updateNode.LEN
	var propRat float64
	updateNode.LEN, propRat = singleBrlenMultiplierProp(updateNode.LEN, chain.STEPLEN)
	llstar := chain.TREELL.CalcCluster(chain, true, cluster)
	lpstar := chain.BRANCHPRIOR.Calc(chain.NODES)
	alpha := math.Exp(lpstar-chain.BRANCHPRIOR.CUR) * math.Exp(llstar-chain.TREELL.CUR) * propRat
	//fmt.Println(llstar, ll, llstar-ll)
	s1 = rand.NewSource(time.Now().UnixNano())
	r1 = rand.New(s1)
	r := r1.Float64()
	if r < alpha {
		//TODO: need to add attribute for cluster-specific acceptance probs in LL struct
		chain.TREELL.CUR = llstar
		chain.BRANCHPRIOR.CUR = lpstar
	} else {
		updateNode.LEN = soldL
	}
}

//this move prunes and regrafts a fossil, creating random variables for the altered branch lengths
//the procedure is basically the same as the SPR move as described in the Yang (2014) Mol. Evol. book
func (chain *MCMC) fossilPlacementUpdate() {
	fn := drawRandomNode(chain.FOS) //draw a random fossil tip
	reattach := drawRandomReattachment(fn, chain.NODES[1:])
	x, p, lastn := PruneFossilTip(fn)
	r := GraftFossilTip(fn.PAR, reattach)
	propRat := r / (x * p)
	//tree.UnmarkAll()
	//lastn.UnmarkToRoot(tree)
	//fn.UnmarkToRoot(tree)
	//reattach.UnmarkToRoot(tree)
	//llstar := WeightedUnrootedLogLike(tree, true, weights)
	llstar := chain.TREELL.Calc(chain.TREE, true)
	//fmt.Println(llstar, chain.TREELL.CUR)
	//llstar1 := WeightedUnrootedLogLike(tree, true, weights)
	//MarkAll(nodes)
	//fmt.Println(llstar, llstar1)
	lpstar := chain.BRANCHPRIOR.Calc(chain.NODES)
	alpha := math.Exp(lpstar-chain.BRANCHPRIOR.CUR) * math.Exp(llstar-chain.TREELL.CUR) * propRat
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r = r1.Float64()
	if r < alpha {
		chain.TREELL.CUR = llstar
		chain.BRANCHPRIOR.CUR = lpstar
	} else { //move fossil back to its previous position and restore old branch lengths
		PruneFossilTip(fn)
		GraftFossilTip(fn.PAR, lastn)
		lastn.LEN = x
		fn.PAR.LEN = p
		//fn.UnmarkToRoot(tree)
		//reattach.UnmarkToRoot(tree)
		//tree.UnmarkAll()
	}
}

func (chain *MCMC) singleBranchLengthUpdate() {
	updateNode := RandomNode(chain.NODES)
	soldL := updateNode.LEN
	//updateNode.UnmarkToRoot(tree)
	var propRat float64
	updateNode.LEN, propRat = singleBrlenMultiplierProp(updateNode.LEN, chain.STEPLEN)
	//updateNode.UnmarkToRoot(tree)
	llstar := chain.TREELL.Calc(chain.TREE, true)
	//MarkAll(nodes)

	lpstar := chain.BRANCHPRIOR.Calc(chain.NODES)

	alpha := math.Exp(lpstar-chain.BRANCHPRIOR.CUR) * math.Exp(llstar-chain.TREELL.CUR) * propRat
	//fmt.Println(llstar, ll, llstar-ll)
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r := r1.Float64()
	if r < alpha {
		chain.TREELL.CUR = llstar
		chain.BRANCHPRIOR.CUR = lpstar
	} else {
		updateNode.LEN = soldL
	}
}

func (chain *MCMC) cladeBranchLengthUpdate() {
	updateNode := RandomInternalNode(chain.NODES)
	updateClade := updateNode.PostorderArray()
	var oldlens []float64
	for _, node := range updateClade {
		oldlens = append(oldlens, node.LEN)
	}
	newlens, propRat := cladeBrlenMultiplierProp(oldlens, chain.STEPLEN)
	for i, node := range updateClade {
		node.LEN = newlens[i]
		node.MRK = false
	}
	//updateNode.UnmarkToRoot(tree)
	llstar := chain.TREELL.Calc(chain.TREE, true)
	lpstar := chain.BRANCHPRIOR.Calc(chain.NODES)

	alpha := math.Exp(lpstar-chain.BRANCHPRIOR.CUR) * math.Exp(llstar-chain.TREELL.CUR) * propRat
	//fmt.Println(llstar, ll, llstar-ll)
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r := r1.Float64()
	if r < alpha {
		chain.TREELL.CUR = llstar
		chain.TREELL.LAST = lpstar
	} else {
		//updateNode.UnmarkToRoot(tree)
		for i, node := range updateClade {
			node.LEN = oldlens[i]
		}
	}
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
