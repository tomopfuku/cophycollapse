package cophymaru

import (
	"gonum.org/v1/gonum/mat"
)

//UVNLike will store all data and params for each cluster
type UVNLike struct {
	Sites []int
	//MeanData *mat.Dense
	DataSum float64
	ProdSum float64 // gives the sum of A*A^T across all sites A in the cluster
	KN      float64
	MuN     float64
	AlphaN  float64
	BetaN   float64
}

//ClusterSampleMean will calculate the sample mean of all the sites in a particular cluster
func (uvn *UVNLike) ClusterSampleMean(aln *CharAlignment) *mat.Dense {
	mean := mat.NewDense(aln.Dim, 1, nil)
	for k := range uvn.Sites {
		mean.Add(mean, aln.MatSites[k])
	}
	recipDim := 1. / float64(len(uvn.Sites))
	mean.Scale(recipDim, mean)
	return mean
}
