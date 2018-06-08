package cophycollapse

import (
	"gonum.org/v1/gonum/mat"
)

//MVNLike will store all data and params for each cluster
type MVNLike struct {
	Sites []int
	//MeanData *mat.Dense
	DataSum *mat.Dense
	ProdSum *mat.Dense // gives the sum of A*A^T across all sites A in the cluster
	kN      float64
	vN      float64
	muN     *mat.Dense
	sN      *mat.Dense
	sNdet   float64
}

//ClusterSampleMean will calculate the sample mean of all the sites in a particular cluster
func (mvn *MVNLike) ClusterSampleMean(aln *CharAlignment) *mat.Dense {
	mean := mat.NewDense(aln.Dim, 1, nil)
	for k := range mvn.Sites {
		mean.Add(mean, aln.MatSites[k])
	}
	recipDim := 1. / float64(len(mvn.Sites))
	mean.Scale(recipDim, mean)
	return mean
}
