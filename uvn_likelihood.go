package cophymaru

//UVNLike will store all data and params for each cluster
type UVNLike struct {
	Sites []int
	//MeanData *mat.Dense
	DataSum []float64
	KN      []float64
	MuN     []float64
	AlphaN  []float64
	BetaN   []float64
}

/*
//ClusterSampleMean will calculate the sample mean of all the sites in a particular cluster
func (uvn *UVNLike) ClusterSampleMean(dm *DistanceMatrix, pairlab int) []float64 {
	dsum := 0.
	for site := range uvn.Sites {
		dsum += dm.MatSites[pairlab][site]
	}
}
*/
