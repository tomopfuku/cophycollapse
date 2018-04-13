package cophymaru

import (
	"fmt"
	"math"
	"math/rand"
	"os"

	"gonum.org/v1/gonum/mat"
)

//DPPGibbs is a struct to store the attributes of a collapsed DPP gibbs sampler for the BM phylo gene expression data
type DPPGibbs struct {
	Aln             *CharAlignment
	Clusters        map[int]*MVNLike
	SiteAssignments map[int]int
	Prior           *GIWPrior
	Gen             int
	PrintFreq       int
	WriteFreq       int
	Threads         int
	Workers         int
	Alpha           float64
}

//InitGibbs will initialize the parameters for the collapsed gibbs sampler.
func InitGibbs(traits map[string][]float64, prior *GIWPrior, gen, print, write, threads, workers int, alpha float64) *DPPGibbs {
	gibbs := new(DPPGibbs)
	aln := InitTraitMatrices(traits)
	gibbs.Aln = aln
	gibbs.Prior = prior
	gibbs.Gen = gen
	gibbs.Alpha = alpha
	gibbs.PrintFreq = print
	gibbs.WriteFreq = write
	gibbs.Threads = threads
	gibbs.Workers = workers
	gibbs.startingClusters()
	return gibbs
}

//Run will run MCMC simulations
func (chain *DPPGibbs) Run() {
	for gen := 0; gen <= chain.Gen; gen++ {
		chain.collapsedGibbsClusterUpdate()
	}
}

func (chain *DPPGibbs) collapsedGibbsClusterUpdate() {
	for k, v := range chain.SiteAssignments {
		chain.collapsedSiteClusterUpdate(k, v)
	}
}

//this will 'reseat' site i
func (chain *DPPGibbs) collapsedSiteClusterUpdate(site int, siteClusterLab int) {
	siteCluster := chain.Clusters[siteClusterLab]
	if len(siteCluster.Sites) == 1 { // delete the cluster containing the current site if it is a singleton
		delete(chain.Clusters, siteClusterLab)
	} else { // delete current site from its cluster
		var newSlice []int
		for _, s := range siteCluster.Sites {
			if s != site {
				newSlice = append(newSlice, s)
			}
		}
		siteCluster.Sites = newSlice
	}
	/*
		for k, v := range chain.Clusters {
			p, _, _, _ := chain.calcClusterProb(k, v, site)
			fmt.Println(p)
		}
	*/
	sum := 0.
	clustProbs := make(map[int]float64)
	clustSn := make(map[int]*mat.Dense)
	clustS := make(map[int]*mat.Dense)
	clustDet := make(map[int]float64)
	for k, v := range chain.Clusters {
		prob, Sn, S, SnD := chain.calcClusterProb(k, v, site)
		sum += prob
		clustProbs[k] = prob
		clustSn[k] = Sn
		clustS[k] = S
		clustDet[k] = SnD
	}
	newClustLab := MaxClustLab(clustProbs) + 1
	newClustProb := chain.Alpha / (float64(chain.Aln.NSites) + chain.Alpha - 1) * chain.Prior.PPDensity
	fmt.Println(chain.Prior.PPDensity)
	clustProbs[newClustLab] = newClustProb
	sum += newClustProb
	r := rand.Float64()
	cumprob := 0.
	newCluster := -1
	for k, v := range clustProbs {
		//clustProbs[k] = v / sum
		cumprob += v / sum
		fmt.Println(cumprob)
		if cumprob > r {
			newCluster = k
			break
		}
	}
	fmt.Println(newCluster)
	os.Exit(0)
	if newCluster < 0 {
		fmt.Println("there was an error assigning a new cluster for site", site)
		os.Exit(0)
	}
	new := chain.Clusters[newCluster]
	new.Sites = append(new.Sites, site)
	new.sN = clustSn[newCluster]
	new.sNdet = clustDet[newCluster]
	new.ProdSum = clustS[newCluster]
	fmt.Println(siteClusterLab, newCluster)
}

// calculates the probability of placing site in clust k
func (chain *DPPGibbs) calcClusterProb(k int, clust *MVNLike, site int) (float64, *mat.Dense, *mat.Dense, float64) {
	N := float64(len(clust.Sites) + 1.)
	kn := (chain.Prior.K0 + N)
	vn := chain.Prior.V0 + N
	knMinusI := (chain.Prior.K0 + (N - 1.))
	vnMinusI := chain.Prior.V0 + (N - 1.)
	DPPprior := (N - 1) / (float64(chain.Aln.NSites) + chain.Alpha - 1.)
	ntax := chain.Aln.Dim
	Sn, S := chain.CalcSn(clust, site, N)
	ppDimExp := float64(-ntax) / 2.
	pp1 := math.Pow(3.141592653589793, ppDimExp)
	pp2 := math.Pow((kn / knMinusI), ppDimExp)
	ppVnExp := vn / 2.
	ppVnMinusIExp := vnMinusI / 2.
	bot := math.Pow(clust.sNdet, ppVnMinusIExp)
	ld, s := mat.LogDet(Sn)
	ld = ld * s
	ld = math.Exp(ld)
	top := math.Pow(ld, ppVnExp)
	pp3 := top / bot
	prob := pp1 * pp2 * pp3
	//fmt.Println(prob, pp1, pp2, pp3, ld)
	t, ts := math.Lgamma(vn / 2.)
	b, bs := math.Lgamma((vn - float64(ntax)) / 2.)
	t = t * float64(ts)
	b = b * float64(bs)
	pp4 := math.Exp(t - b)
	prob = prob * pp4
	prob = prob * DPPprior
	//fmt.Println(prob, pp1, pp2, pp3, pp4)
	return prob, Sn, S, ld
}

func calcGIWPostParams() {
}

//CalcMn will calculate the mean vector of the posterior for a single cluster after a new site has been added
func (chain *DPPGibbs) CalcMn(clust *MVNLike, N, kn float64, newSiteData *mat.Dense) *mat.Dense {
	//clust := chain.Clusters[k]
	r, _ := clust.DataSum.Dims()
	avgX := mat.NewDense(r, 1, nil)
	avgX.Add(clust.DataSum, newSiteData)
	avgX.Scale((1. / N), avgX)
	top := mat.NewDense(r, 1, nil)
	top.Scale(N, avgX)
	top.Add(chain.Prior.Mu0K0, top)
	top.Scale(kn, top)
	//matPrint(top)
	return top
}

//CalcSn will calculate the VCV of the posterior distribution
func (chain *DPPGibbs) CalcSn(clust *MVNLike, site int, N float64) (*mat.Dense, *mat.Dense) {
	ntax := chain.Aln.Dim
	newSiteProduct := chain.Aln.MatSitesProd[site]
	S := mat.NewDense(ntax, ntax, nil)
	S.Add(clust.ProdSum, newSiteProduct)
	S0 := chain.Prior.S0
	Sn := mat.NewDense(ntax, ntax, nil)
	Sn.Add(S0, S)
	mu0Mult := chain.Prior.PostPredProd
	Sn.Add(Sn, mu0Mult)
	kN := clust.kN + 1.
	divKn := 1. / kN
	muN := chain.CalcMn(clust, N, kN, chain.Aln.MatSites[site])
	muNMult := mat.NewDense(ntax, ntax, nil)
	muNMult.Product(muN, muN.T())
	muNMult.Scale(divKn, muNMult)
	Sn.Sub(Sn, muNMult)
	return Sn, S
}

func (chain *DPPGibbs) startingClusters() {
	clus := make(map[int]*MVNLike)
	lab := 0
	siteClust := make(map[int]int)
	for k := range chain.Aln.MatSites {
		curMVN := new(MVNLike)
		curMVN.Sites = []int{k}
		curMVN.DataSum = chain.Aln.MatSites[k]
		curMVN.ProdSum = chain.Aln.MatSitesProd[k]
		curMVN.kN = chain.Prior.K0 + 1.
		curMVN.vN = chain.Prior.V0 + 1.
		chain.initClusterMatrixParams(curMVN)
		clus[lab] = curMVN
		siteClust[k] = lab
		lab++
	}
	chain.Clusters = clus
	chain.SiteAssignments = siteClust
}

func (chain *DPPGibbs) initClusterMatrixParams(clust *MVNLike) {
	r, _ := clust.DataSum.Dims()
	muN := mat.NewDense(r, 1, nil)
	priorMult := chain.Prior.Mu0K0
	top := clust.DataSum
	top.Add(priorMult, top)
	muN.Scale((1. / clust.kN), top)
	//muN is good
	sN := mat.NewDense(r, r, nil)
	sN.Add(chain.Prior.S0, clust.ProdSum)
	sN.Add(sN, chain.Prior.PostPredProd)
	muNprod := mat.NewDense(r, r, nil)
	muNprod.Product(muN, muN.T())
	muNprod.Scale((1. / clust.kN), muNprod)
	sN.Sub(sN, muNprod)
	clust.sN = sN
	clust.muN = muN
	ld, s := mat.LogDet(sN)
	ld = ld * s
	ld = math.Exp(ld)
	clust.sNdet = ld
}
