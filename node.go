package cophymaru

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
)

//Node is a node struct to represent rooted and unrooted phylogenetic trees
type Node struct {
	PAR    *Node
	CHLD   []*Node
	NAME   string
	LEN    float64
	PRNLEN float64
	CONTRT []float64
	//VISITED bool
	MRK         bool
	MIS         []bool //this is a slice of bools that indicates whether the index is missing from CONTRT
	LL          []float64
	CONPRNLEN   []float64
	ClustLEN    map[int]float64
	ClustPrnLen [][]float64
}

//PostorderArray will return an array of all the nodes in the tree in Postorder
func (n *Node) PostorderArray() (ret []*Node) {
	var buffer []*Node
	for _, cn := range n.CHLD {
		for _, cret := range cn.PreorderArray() {
			buffer = append(buffer, cret)
		}
	}
	buffer = append(buffer, n)
	ret = buffer
	return
}

//PreorderArray will return a preordered array of all the nodes in a tree
func (n *Node) PreorderArray() (ret []*Node) {
	var buffer []*Node
	buffer = append(buffer, n)
	for _, cn := range n.CHLD {
		for _, cret := range cn.PreorderArray() {
			buffer = append(buffer, cret)
		}
	}
	ret = buffer
	return
}

//Newick will return a newick string representation of the tree
func (n Node) Newick(bl bool) (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.CHLD {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Newick(bl))
		if bl == true {
			s := strconv.FormatFloat(cn.LEN, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.CHLD)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.NAME)
	ret = buffer.String()
	return
}

//AddChild will add a child to a node
func (n *Node) AddChild(c *Node) error {
	for _, ch := range n.CHLD {
		if c == ch {
			fmt.Println("You are trying to reroot on the current node!")
		}
	}
	n.CHLD = append(n.CHLD, c)
	c.PAR = n
	return nil
}

//RemoveChild will remove a chidl from the slice of children associated with a node
func (n *Node) RemoveChild(c *Node) {
	var newCHLD []*Node
	for _, ch := range n.CHLD {
		if ch != c {
			newCHLD = append(newCHLD, ch)
		}
		n.CHLD = nil
		n.CHLD = newCHLD
		c.PAR = nil
	}
}

//NNodes is a helper method that will return the number of internal nodes descending from n (including n)
func (n *Node) NNodes(count *int) {
	*count++
	for _, ch := range n.CHLD {
		ch.NNodes(count)
	}
}

//Reroot reroots all the nodes represented in a graph on n
func (n *Node) Reroot(oldroot *Node) *Node {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	nnodes := 0
	oldroot.NNodes(&nnodes)
	var pathnodes = make([]*Node, nnodes)
	//var pathnodes []*Node
	curnode := n
	pathlen := 0 //this will count the number of nodes between the newroot and the oldroot
	for ind := range pathnodes {
		pathnodes[ind] = curnode
		//pathnodes = append(pathnodes,curnode)
		if curnode == oldroot {
			break
		}
		pathlen++
		curnode = curnode.PAR
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.RemoveChild(newpar)
		newpar.AddChild(curnode)
		curnode.LEN = newpar.LEN
	}
	//curnode = nil
	//newpar = nil
	n.LEN = 0.0
	return n
}

//UnmarkToRoot will mark all of the nodes betwen the current node and the root to be recalculated
func (n *Node) UnmarkToRoot(oldroot *Node) {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	curnode := n
	pathlen := int(0)
	for {
		pathlen++
		curnode = curnode.PAR
		curnode.MRK = false
		if curnode == oldroot {
			break
		}
		if pathlen > 2000000 {
			fmt.Println("something is wrong with the tree. could not walk back to root.")
			os.Exit(0)
		}
	}
}

//UnmarkAll will unmark all of the nodes on a tree
func (n *Node) UnmarkAll() {
	nodes := n.PreorderArray()
	for _, node := range nodes {
		node.MRK = false
	}
}

//MarkAll will unmark all of the nodes on a tree
func (n *Node) MarkAll() {
	nodes := n.PreorderArray()
	for _, node := range nodes {
		node.MRK = true
	}
}
