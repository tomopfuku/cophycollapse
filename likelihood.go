package cophymaru

//LL is a struct for likelihood functions/calculations
type LL struct {
	MULTI       bool
	WORKERS     int
	SITEWEIGHTS []float64
	CUR         float64
	LAST        float64
}

//CalcCluster will calculate the log-likelihood of a single cluster
func (ll *LL) CalcCluster(chain *MCMC, startFresh bool, cluster int) float64 {
	return ClusterLogLike(chain, cluster, startFresh, ll.WORKERS)
}

//Calc will calculate the log-likelihood
func (ll *LL) Calc(tree *Node, startFresh bool) float64 {
	if ll.MULTI == true {
		return WeightedUnrootedLogLikeParallel(tree, startFresh, ll.SITEWEIGHTS, ll.WORKERS)
	}
	return WeightedUnrootedLogLike(tree, startFresh, ll.SITEWEIGHTS)
}

//InitLL will initialize the likelihood struct
func InitLL(multithread bool, workers int, weights []float64) *LL {
	ll := new(LL)
	if multithread == true {
		ll.MULTI = true
		ll.WORKERS = workers
		ll.SITEWEIGHTS = weights
	} else {
		ll.MULTI = false
		ll.SITEWEIGHTS = weights
	}
	return ll
}
