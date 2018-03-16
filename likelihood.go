package cophymaru

//LL is a struct for likelihood functions/calculations
type LL struct {
	MULTI       bool
	WORKERS     int
	SITEWEIGHTS []float64
	CUR         float64
	LAST        float64
	CLUSTCUR    map[int]float64
	CLUSTLAST   map[int]float64
}

//CalcCluster will calculate the log-likelihood of a single cluster
func (ll *LL) CalcCluster(chain *MCMC, startFresh bool, cluster int) float64 {
	nsites := len(chain.CLUSTERSET[cluster])
	if nsites != 1 {
		return ClusterLogLikeParallel(chain, cluster, startFresh, ll.WORKERS)
	}
	return ClusterLogLike(chain, cluster, startFresh)
}

//Calc will calculate the log-likelihood
func (ll *LL) Calc(tree *Node, startFresh bool) float64 {
	if ll.MULTI == true {
		return WeightedUnrootedLogLikeParallel(tree, startFresh, ll.SITEWEIGHTS, ll.WORKERS)
	}
	return WeightedUnrootedLogLike(tree, startFresh, ll.SITEWEIGHTS)
}

//CalcCombinedClusterLL will calculate the total tree log-likelihood
func (ll *LL) CalcCombinedClusterLL(chain *MCMC) (treelik float64) {
	treelik = 0.
	for c := range chain.UNIQUEK {
		AssignClustLens(chain, c)
		cur := ll.CalcCluster(chain, true, c)
		treelik += cur
	}
	return
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
