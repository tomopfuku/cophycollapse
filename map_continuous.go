package cophycollapse

import (
	"fmt"
	//"io/ioutil"
	"strconv"
	"strings"
)

//this will be used to extract the phylip header so that we can allocate according to ntax and ntraits
func readPhylipHeader(traitfl string) (traits []string, ntax int, ntraits int) {
	lines := ReadLine(traitfl)
	//var ntax,ntraits int
	var err error
	ss := strings.Split(lines[0], "\t")
	ntax, err = strconv.Atoi(ss[0])
	if err != nil {
		fmt.Println("had trouble reading in the PHYLIP header. Check to make sure that it is correct")
		fmt.Println(lines[0])
	}
	ntraits, err = strconv.Atoi(ss[1])
	if err != nil {
		fmt.Println("had trouble reading in the PHYLIP header. Check to make sure that it is correct.")
		fmt.Println(lines[0])
	}
	traits = lines[1:]
	return
}

// ReadContinuous will read in a PHYLIP-formatted trait file and return a map
func ReadContinuous(traitfl string) (map[string][]float64, int, int) {
	traitlines, ntax, ntraits := readPhylipHeader(traitfl)
	tm := make(map[string][]float64, ntax)
	//lnum := 0
	//var ntax,ntraits int
	var err error
	var curtr float64
	for _, line := range traitlines {
		if line == "" {
			continue
		}
		ss := strings.Split(line, "\t")
		curtax := ss[0]
		var curtraits []float64
		if len(ss[1:]) != ntraits {
			fmt.Println("the number of traits specified in the header differs from the number present in this line:")
			fmt.Println(line)
		}
		for _, v := range ss[1:] {
			if v != "?" {
				curtr, err = strconv.ParseFloat(v, 64)
				if err != nil {
					fmt.Println("couldn't read in the trait file starting at this line:")
					fmt.Println(line)
				}
			} else {
				curtr = -1000000.0
			}
			curtraits = append(curtraits, curtr)
		}
		tm[curtax] = curtraits
	}
	return tm, ntax, ntraits
}

//MakeMissingDataSlice will intialize the MIS attribute of node t by identifying traits with LARGE (ie. missing) values and updating MIS accordingly.
func MakeMissingDataSlice(t *Node) {
	for _, tr := range t.CONTRT { //create slice marking missing data
		if tr == -1000000.0 {
			t.MIS = append(t.MIS, true)
		} else {
			t.MIS = append(t.MIS, false)
		}
	}
}

//MapContinuous maps the traits contained with in a (golang) map to the tips of a tree and initializes slices of the same length for the internal nodes
func MapContinuous(t *Node, traits map[string][]float64, ntraits int) {
	var z float64
	z = 0.000000000
	if len(t.CHLD) == 0 {
		t.CONTRT = traits[t.NAME]
		t.LL = append(t.LL, z)
		MakeMissingDataSlice(t)
	} else {
		count := 0
		for {
			if count == ntraits {
				break
			}
			t.CONTRT = append(t.CONTRT, z)
			t.MIS = append(t.MIS, false)
			t.LL = append(t.LL, z)
			count++
		}
	}
	for _, chld := range t.CHLD {
		MapContinuous(chld, traits, ntraits)
	}
}
