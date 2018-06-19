package cophycollapse

import (
	"gonum.org/v1/gonum/stat/distuv"
)

//InitNGPrior will initialize the prior on the  matrix and BM mean
func InitNGPrior(mu0, k0, alpha0, beta0 float64) *NormalGammaPrior {
	p := new(NormalGammaPrior)
	p.Mu0 = mu0
	p.K0 = k0
	p.Alpha0 = alpha0
	p.Beta0 = beta0
	stt := new(distuv.StudentsT)
	stt.Mu = p.Mu0
	stt.Sigma = (p.Beta0 * (p.K0 + 1.)) / (p.Alpha0 * p.K0)
	stt.Nu = 2 * p.Alpha0
	p.PPDensity = stt
	return p
}

//NormalGammaPrior is an object for the Gausian inverse Wishart prior used as conjugate prior for brownain phylogenetic VCV likelihood
type NormalGammaPrior struct {
	Mu0    float64 // hyperparameter for expected mean pairwise distance
	K0     float64 // hyperparameter on the precision (1/variance) of the prior mean
	Alpha0 float64 // shape hyperpameter on the variance of the distance estimates
	Beta0  float64 // scale hyperparemeter on the variance of the distance estimates
	//PostPredProd float64 // this will store Mu0*Mu0^T*K0 (for use in calculating the posterior parameters)
	//Mu0K0        *mat.Dense
	PPDensity *distuv.StudentsT
}

// NGPost stores values for the posterior distribution calculated analytically from GIWPrior and the phylo BM likelihood
type NGPost struct {
	MuN    float64 // parameters are same as prior
	KN     float64
	AlphaN float64
	BetaN  float64
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

/*/GIWStartingSampleMean will calculate the starting mean value as just the raw sample mean from all traits
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
*/
