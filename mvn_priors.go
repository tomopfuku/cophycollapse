package cophycollapse

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

//InitVCVPrior will initialize the prior on the VCV matrix and BM mean
func InitVCVPrior(mu0 *mat.Dense, k0 float64, s0 *mat.Dense, v0 float64) *GIWPrior {
	prior := new(GIWPrior)
	prior.Mu0 = mu0
	prior.K0 = k0
	prior.S0 = s0
	prior.V0 = v0
	dim, _ := mu0.Dims()
	prod := mat.NewDense(dim, dim, nil)
	prod.Product(mu0, mu0.T())
	prod.Scale(k0, prod)
	prior.PostPredProd = prod
	m0k0 := mat.NewDense(dim, 1, nil)
	m0k0.Scale(k0, mu0)
	prior.Mu0K0 = m0k0
	D := float64(dim) / 2.
	s0det, s := mat.LogDet(s0)
	s0det = math.Exp(s0det * s)
	fmt.Println(math.Pow(s0det, -0.5))
	pp := math.Pow(3.141592653589793, -D) * (math.Pow(k0, D) / math.Pow((k0+1.), D)) * math.Pow(s0det, -0.5)
	rprod := 0.
	for i := 0; i <= dim; i++ {
		//rprod = rprod * (math.Gamma((v0+2-float64(i))/2.) / math.Gamma((v0+1-float64(i))/2.))
		top, s := math.Lgamma((v0 + 2 - float64(i)) / 2.)
		top = top * float64(s)
		bot, s := math.Lgamma((v0 + 1 - float64(i)) / 2.)
		bot = bot * float64(s)
		rprod = rprod + (bot - bot)
		//fmt.Println(i, rprod)
	}
	pp = pp * math.Exp(rprod)
	prior.PPDensity = pp
	return prior
}

//GIWPrior is an object for the Gausian inverse Wishart prior used as conjugate prior for brownain phylogenetic VCV likelihood
type GIWPrior struct {
	Mu0          *mat.Dense // prior on multivariate normal mean (ancestral state at root of tree)
	K0           float64    // belief in Mu0
	S0           *mat.Dense // prior on multivariate normal VCV
	V0           float64    // belief in S0
	PostPredProd *mat.Dense // this will store Mu0*Mu0^T*K0 (for use in calculating the posterior parameters)
	Mu0K0        *mat.Dense
	PPDensity    float64
}

// GIWPost stores values for the posterior distribution calculated analytically from GIWPrior and the phylo BM likelihood
type GIWPost struct {
	MuN *mat.Dense // parameters are same as prior
	KN  float64
	S   *mat.Dense
	VN  float64
}

/*
func (p *GIWPost) calcPostParams(prior *GIWPrior, like *MVNLike) {
	Dim := float64(like.Dim)
	kN := Dim + prior.K0
	vN := Dim + prior.V0
	p.KN = kN
	p.VN = vN
}
*/

//GIWStartingSampleMean will calculate the starting mean value as just the raw sample mean from all traits
func GIWStartingSampleMean(traits map[string][]float64) (mu0 *mat.Dense) {
	count := 0.
	nobs := 0.
	for k := range traits {
		for _, t := range traits[k] {
			count += t
			nobs += 1.
		}
	}
	avg := count / nobs
	dim := len(traits)
	mu0 = mat.NewDense(dim, 1, nil)
	for i := 0; i < dim; i++ {
		mu0.Set(i, 0, avg)
	}
	return
}
