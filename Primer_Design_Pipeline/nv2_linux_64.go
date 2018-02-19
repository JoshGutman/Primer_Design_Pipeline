package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"regexp"
	"strings"
	"sync"
)

type Genome struct {
	name string
	seq  string
}

type Match struct {
	genome *Genome
	index  int
}

type Amplicon struct {
	genome  *Genome
	forward int
	reverse int
	seq     string
}

func main() {

	g := flag.String("g", "", "[REQUIRED] Path to .fasta file")
	f := flag.String("f", "", "[REQUIRED] Forward sequence")
	r := flag.String("r", "", "[REQUIRED] Reverse sequence")
	m := flag.Int("m", (1<<31)-1, "Max length of amplicon")

	flag.Parse()

	if *g == "" || *f == "" || *r == "" {
		fmt.Println("Incorrect usage. Program exiting...")
		os.Exit(1)
	}

	fasta := readFasta(*g)
	genomes := makeGenomes(fasta)
	maxLen := *m
	
	forwards := expandDegenerateSequence([]byte(*f))
	reverses := expandDegenerateSequence([]byte(reverseComplement(*r)))
	antiForwards := expandDegenerateSequence([]byte(reverseComplement(*f)))
	antiReverses := expandDegenerateSequence([]byte(*r))
	

	fChan := make(chan []Match, (len(forwards)+1)*len(genomes))
	rChan := make(chan []Match, (len(reverses)+1)*len(genomes))
	afChan := make(chan []Match, (len(antiForwards)+1)*len(genomes))
	arChan := make(chan []Match, (len(antiReverses)+1)*len(genomes))

	chanArray := []chan []Match{fChan, rChan, afChan, arChan}
	
	var wg sync.WaitGroup

	
	for _, genome := range genomes {
	
	
		
		for _, f := range forwards {
			wg.Add(1)
			go func(seq string, g Genome) {
				defer wg.Done()
				fChan <- findMatches(regexp.MustCompile(seq), &g)
			}(string(f), genome)
		}

		for _, r := range reverses {
			wg.Add(1)
			go func(seq string, g Genome) {
				defer wg.Done()
				rChan <- findMatches(regexp.MustCompile(seq), &g)
			}(string(r), genome)
		}
		
		for _, af := range antiForwards {
			wg.Add(1)
			go func(seq string, g Genome) {
				defer wg.Done()
				afChan <- findMatches(regexp.MustCompile(seq), &g)
			}(string(af), genome)
		}
		
		for _, ar := range antiReverses {
			wg.Add(1)
			go func(seq string, g Genome) {
				defer wg.Done()
				arChan <- findMatches(regexp.MustCompile(seq), &g)
			}(string(ar), genome)
		}
		
	}

	wg.Wait()

	for i, _ := range chanArray {
		close(chanArray[i])
	}
	
	
	fMatches := make(map[string][]Match)
	rMatches := make(map[string][]Match)
	afMatches := make(map[string][]Match)
	arMatches := make(map[string][]Match)
	
	matchArray := []map[string][]Match{fMatches, rMatches, afMatches, arMatches}
	
	for i, _ := range chanArray {
		for subArray := range chanArray[i] {
			for _, match := range subArray {
				matchArray[i][match.genome.name] = append(matchArray[i][match.genome.name], match)
			}
		}
	}

	
	amps := make(map[string][]Amplicon)
	antiAmps := make(map[string][]Amplicon)
	var min int
	
	if len(fMatches) < len(rMatches) {
		min = 0
	} else {
		min = 1
	}
	for k, _ := range matchArray[min] {
		amps[k] = findAmpliconsInRange(fMatches[k], rMatches[k], maxLen, len(*r))
	}
	
	if len(afMatches) < len(arMatches) {
		min = 2
	} else {
		min = 3
	}
	for k, _ := range matchArray[min] {
		antiAmps[k] = findAmpliconsInRange(arMatches[k], afMatches[k], maxLen, len(*r))
	}
	
	for _, ampArray := range amps {
		for i := 0; i < len(ampArray); i++ {
			ampArray[i].getAmpliconSeq(len(*r))
		}
	}

	for _, ampArray := range antiAmps {
		for i := 0; i < len(ampArray); i++ {
			ampArray[i].getAmpliconSeq(len(*r))
		}
	}
	
	for _, ampArray := range amps {
		for _, amp := range ampArray {
			fmt.Printf("%d\t%s\t%s\n", len(amp.seq), amp.seq, amp.genome.name)
		}
	}
	
	for _, ampArray := range antiAmps {
		for _, amp := range ampArray {
			fmt.Printf("%d\t%s\t%s\n", len(amp.seq), amp.seq, amp.genome.name)
		}
	}
}

func (a *Amplicon) getAmpliconSeq(reverseLen int) {
	a.seq = (a.genome.seq)[a.forward:(a.reverse + reverseLen)]
}

func findAmpliconsInRange(forwardMatches, reverseMatches []Match, size, reverseLen int) []Amplicon {
	var out []Amplicon

	for _, fm := range forwardMatches {
		for _, rm := range reverseMatches {
			ampLen := (rm.index - fm.index + reverseLen)
			if ampLen <= size && ampLen > 0 {
				out = append(out, Amplicon{fm.genome, fm.index, rm.index, ""})
			}
		}
	}

	return out
}


func findMatches(re *regexp.Regexp, g *Genome) []Match {
	var out []Match
	matches := re.FindAllStringIndex((*g).seq, -1)
	for _, m := range matches {
		out = append(out, Match{g, m[0]})
	}
	return out
}

func readFasta(infile string) [][]byte {

	// 1. Read in entire file
	// 2. Split file based on ">"
	// 3. For each fasta, split based on "\n" only one time
	// 4. Ignore element 0 (will be genome name. Element 1 will be sequence)
	// TODO:
	// 5. For each sequence, do main() (Will have to encapsulate some work in separate function) (Threaded?)

	contents, err := ioutil.ReadFile(infile)
	if err != nil {
		log.Fatal(err)
	}

	g := bytes.Split(contents, []byte{'>'})

	return g[1:]
}

func makeGenomes(gnms [][]byte) []Genome {
	var genomes []Genome
	for _, g := range gnms {
		info := bytes.SplitN(g, []byte{'\n'}, 2)
		genomes = append(genomes, Genome{string(bytes.Replace(info[0], []byte{'\n'}, []byte{}, -1)), string(bytes.Replace(info[1], []byte{'\n'}, []byte{}, -1))})
	}
	return genomes
}

func reverseComplement(seq string) string {
	ambiguityCodes := map[uint8]uint8{
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'M': 'K',
		'R': 'Y',
		'W': 'W',
		'S': 'S',
		'Y': 'R',
		'K': 'M',
		'V': 'B',
		'H': 'D',
		'D': 'H',
		'B': 'V',
		'N': 'N',
		'.': '.',
		'-': '-',
		'^': '^'}

	seq = strings.ToUpper(seq)
	out := ""
	for i := len(seq) - 1; i >= 0; i-- {
		out += string(ambiguityCodes[seq[i]])
	}
	return out
}

func expandDegeneratePosition(primers [][]byte, position int, l ...byte) [][]byte {
	for j := range primers {
		primers[j][position] = l[0]
	}
	primers_len := len(primers)
	for _, m := range l[1:] {
		for j := 0; j < primers_len; j++ {
			p := make([]byte, len(primers[j]))
			copy(p, primers[j])
			p[position] = m
			primers = append(primers, p)
		}
	}
	return primers
}

func expandDegenerateSequence(sequence []byte) [][]byte {
	var primers [][]byte

	if len(sequence) == 0 {
		return nil
	}

	primer := make([]byte, len(sequence))
	copy(primer, sequence)

	primers = append(primers, primer)
	for i, nt := range sequence {
		switch nt {
		default:
			fmt.Println("Error")
		case 'A', 'C', 'G', 'T', 'U':
			primers = expandDegeneratePosition(primers, i, nt)
		case 'W':
			primers = expandDegeneratePosition(primers, i, 'A', 'T')
		case 'S':
			primers = expandDegeneratePosition(primers, i, 'G', 'C')
		case 'M':
			primers = expandDegeneratePosition(primers, i, 'A', 'C')
		case 'K':
			primers = expandDegeneratePosition(primers, i, 'G', 'T')
		case 'R':
			primers = expandDegeneratePosition(primers, i, 'A', 'G')
		case 'Y':
			primers = expandDegeneratePosition(primers, i, 'C', 'T')
		case 'B':
			primers = expandDegeneratePosition(primers, i, 'C', 'G', 'T')
		case 'D':
			primers = expandDegeneratePosition(primers, i, 'A', 'G', 'T')
		case 'H':
			primers = expandDegeneratePosition(primers, i, 'A', 'C', 'T')
		case 'V':
			primers = expandDegeneratePosition(primers, i, 'A', 'C', 'G')
		case 'N':
			primers = expandDegeneratePosition(primers, i, 'G', 'A', 'T', 'C')
		case '.', '-':
			panic(fmt.Errorf("%q is a valid IUPAC Nucleotide Code, but not a valid nucleotide at position %d in sequence %s", nt, i+1, sequence))
		}
	}
	return primers
}
