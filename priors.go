package cophymaru

import (
	"fmt"
	"math"
	"math/big"
	"os"
)

//BranchLengthPrior is a struct for specifying and calculating different priors
type BranchLengthPrior struct {
	TYPE      string
	ALPHA     float64
	BETA      float64
	NTIPS     int
	SGAMMA    float64
	SFACT     float64
	CUR       float64
	LAST      float64
	CLUSTCUR  map[int]float64
	CLUSTLAST map[int]float64
}

//Calc will return the log prior probability of the branch length prior
func (pr *BranchLengthPrior) Calc(nodes []*Node) float64 {
	if pr.TYPE == "1" {
		return ExponentialBranchLengthLogPrior(nodes[1:], pr.BETA)
	} else if pr.TYPE == "2" {
		return DirichletBranchLengthLogPrior(nodes[1:], pr)
	}
	return 0.0
}

func logFactorial(val int) (x float64) {
	x = 0.
	for i := 1; i <= val; i++ {
		x += math.Log(float64(i))
	}
	return
}

//InitializePrior will set up a new instance of a prior
func InitializePrior(priorType string, nodes []*Node) *BranchLengthPrior {
	pr := new(BranchLengthPrior)
	pr.TYPE = priorType
	if pr.TYPE == "0" {
		pr.BETA = 0.
	} else if pr.TYPE == "1" {
		pr.BETA = 1.0
	} else if pr.TYPE == "2" {
		T := TreeLength(nodes)
		pr.ALPHA = 1.0
		//pr.BETA = 10.0
		pr.BETA = T
		ntips := 0
		for _, n := range nodes {
			if len(n.CHLD) == 0 {
				ntips++
			}
		}
		x := (2 * ntips) - 4
		logfact := logFactorial(x)
		pr.NTIPS = ntips
		pr.SFACT = logfact
		pr.SGAMMA = math.Gamma(pr.ALPHA)

	} else {
		fmt.Println("PLEASE SPECIFY VALID OPTION FOR PRIOR. TYPE maru -h TO SEE OPTIONS")
		os.Exit(0)
	}
	return pr
}

func normalPDF(p float64, mean float64, sd float64) float64 {
	prob := (1.0 / math.Sqrt(2.0*math.Pi*sd)) * math.Exp(-(math.Pow(p-mean, 2.) / (2 * sd)))
	return prob
}

func gammaPDF(p, alpha, beta float64) float64 {
	prob := (math.Pow(alpha, beta) / math.Gamma(alpha)) * math.Exp(-beta*p) * math.Pow(p, (alpha-1))
	return prob
}

//GammaTreeLengthPrior will calculate the prior probability of the current tree length
func GammaTreeLengthPrior(nodels []*Node, alpha, beta float64) float64 {
	tl := TreeLength(nodels)
	return gammaPDF(tl, alpha, beta)
}

//DirichletBranchLengthLogPrior will return the log prior probability of all branch lengths in a trejke
func DirichletBranchLengthLogPrior(nodes []*Node, pr *BranchLengthPrior) float64 {
	//NOTE: this is incorrect right now-- need to fix. the big exponent is calculated by rounding to ints at moment-- need to make a big.Float.Pow
	T := TreeLength(nodes)
	x := (2. * float64(pr.NTIPS)) - 4.
	pow := pr.ALPHA - 1. - x
	//prob := math.Log(math.Pow(pr.BETA, pr.ALPHA) / pr.SGAMMA)
	prob := (pr.ALPHA * math.Log(pr.BETA)) - math.Log(pr.SGAMMA)
	prob += (-pr.BETA * T) //+	math.Log(math.Pow(T, (pr.ALPHA - 1. - x))) // tranlated into log math from version in yang book
	prob += (pow * math.Log(T))
	prob = prob + pr.SFACT
	//fmt.Println(logfact, pr.SFACT)
	return prob
}

func factorial(val *big.Int) *big.Int {
	result := new(big.Int)
	switch val.Cmp(&big.Int{}) {
	case -1, 0:
		result.SetInt64(1)
	default:
		result.Set(val)
		var one big.Int
		one.SetInt64(1)
		result.Mul(result, factorial(val.Sub(val, &one)))
	}
	return result
}

//this calculates the probability of drawing parameter p from an exponential distribution with shape parameter == lambda
func exponentialPDF(p float64, beta float64) float64 {
	prob := (1.0 / beta) * math.Exp(-(p / beta))
	return prob
}

//ExponentialBranchLengthLogPrior will give the exponential prior probability of observing a particular set of branch lengths on a tree, given shape lambda
func ExponentialBranchLengthLogPrior(nodels []*Node, beta float64) (prob float64) {
	var tempp float64
	for _, n := range nodels {
		tempp = exponentialPDF(n.LEN, beta)
		prob += math.Log(tempp)
	}
	return
}

//EBExponentialBranchLengthLogPrior will give the exponential probability of each branch length when exponential mean == ML branch length
func EBExponentialBranchLengthLogPrior(nodels []*Node, means []float64) (prob float64) {
	var tempp float64
	for i, n := range nodels {
		tempp = exponentialPDF(n.LEN, means[i])
		prob += math.Log(tempp)
	}
	return
}

//EBNormalBranchLengthLogPrior will calculate the prior probability of a branch length from a normal distribution centered on the ML branch length estimate
func EBNormalBranchLengthLogPrior(nodels []*Node, means []float64) (prob float64) {
	var tempp float64
	for i, n := range nodels {
		tempp = normalPDF(n.LEN, means[i], 0.05)
		prob += math.Log(tempp)
	}
	return
}
