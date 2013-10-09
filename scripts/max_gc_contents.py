#!/usr/bin/env python
import sys

from util import dna_stats

if __name__ == "__main__":
    f = open(sys.argv[1])
    name = ""
    contents = ""
    max_percent = 0
    max_name = ""
    
    all_lines = f.readlines()
    all_lines.append(">")
    for l in all_lines:
        s = l.strip()
        if s.startswith(">"):
            if name != "" and contents != "":
                d = dna_stats.DNAStats(contents, name)
                pct = d.gc_contents()
                if pct > max_percent:
                    max_percent = pct
                    max_name = name
            name = s[1:]
            contents = ""
        else:
            contents += s
    
    print max_name
    print max_percent * 100
