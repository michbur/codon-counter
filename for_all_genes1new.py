from __future__ import division
from collections import OrderedDict

# wczytuje wszystkie sekwencje genow
def readingfile(filename):
        f = open(filename)
        order = []
        sequences = OrderedDict()
        seq = []
        for line in f:
                if line.startswith('>'):
                        name = line[1:].rstrip('\n')
                        name = name.replace('_', ' ')
                        order.append(name)
                        sequences[name] = ''
                else:
                        sequences[name] += line.rstrip('\n').rstrip('*')
        for key in sequences:
                seq.append(sequences[key])
        print "%d znalezionych sekwencji" % len(seq)
        return seq

# dzieli sekwencje na kodony
def SplitOrfsIntoCodons(seqs):
    list_of_split_orfs = []
    for seq in seqs:
        split_ORF = [seq[i:i+3] for i in range(0, len(seq), 3)]
        list_of_split_orfs.append(split_ORF)
    return list_of_split_orfs


DangerousCodons = ["TGG", "TAC", "TAT", "TCA", "TTA", "TGC", "TGT", "GAA", "GAG", "AAA", "AAG", "CAA", "CAG", "TCG", "TTG", "AGA", "CGA", "GGA"]
With_Alt = ["TCA", "TCG", "TTA", "TTG", "AGA", "CGA", "GGA"]
Double_danger = ['TGG', 'TAC', 'TAT', 'TCA', 'TTA']
No_Alt = ['TGG', 'TGT', 'TGC', 'TAC', 'TAT', 'CAA', 'CAG', 'GAA', 'GAG', 'AAA', 'AAG']

# podlicza ilosc wystepowania danych kodonow w 3 czesciach genu
def findandcountcodons(orf, dic):  
        list_of_index = [i for i, j in enumerate(orf) if j in dic]
        first = []
        second = []
        third = []
        a = len(orf)
        b = int((1/3) * a)
        c = int((2/3) * a)
        freq_1 = 0
        freq_2 = 0
        freq_3 = 0
        for index in list_of_index:     
            if index in xrange(0, b):
                first.append(index)
            elif index in xrange(b,c):
                second.append(index)
            else:
                third.append(index)          
            freq_1 = len(first) / len(xrange(0,b))
            freq_2 = len(second) / len(xrange(b,c))
            freq_3 = len(third) / len(xrange(c,a))
        return freq_1, freq_2, freq_3

# wykonuje poprzednia funkcje dla wszystkich genow
def analyzelistoforfs(orf_list, codons):
    results = []
    for eachorf in orf_list:
        s = findandcountcodons(eachorf, codons)
        results.append(s)
    return results

# podlicza czestosc wystepowania kodonow danej grupy w danej czesci genu (dla wszystkich genow)
def countall(orfs):
    part1 = 0
    part2 = 0
    part3 = 0
    a = len(orfs)
    for orf in orfs:
        part1 = part1 + orf[0]
        part2 = part2 + orf[1]
        part3 = part3 + orf[2]
    part1 = part1 / a
    part2 = part2 / a
    part3 = part3 / a
    return part1, part2, part3
        
# wczytuje wszystkie sekwencje genow z pliku
a = readingfile('pr.sek.txt')

# dzieli wszystkie sekwencje na kodony
splitORF = SplitOrfsIntoCodons(a)

# funkcje wykonywane dla wszystkich grup kodonow
b = analyzelistoforfs(splitORF, DangerousCodons)
c = analyzelistoforfs(splitORF, With_Alt)
d = analyzelistoforfs(splitORF, No_Alt)
e = analyzelistoforfs(splitORF, Double_danger)
count1 = countall(b)
count1 = str(count1)

count2= countall(c)
count2 = str(count2)

count3 = countall(d)
count3 = str(count3)

count4 = countall(e)
count4 = str(count4)

# zapisywanie wynikow w plikach 
out = open("beg_mid_end_all_all1.txt", "w")
out.write(count1)
out.close()

out2 = open("beg_mid_end_syn_all1.txt", "w")
out2.write(count2)
out2.close()

out3 = open("beg_mid_end_noalt_all1.txt", "w")
out3.write(count3)
out3.close()


out4 = open("beg_mid_end_double_all1.txt", "w")
out4.write(count4)
out4.close()




# do testow
csv = open('beg_mid_end_all_all1_p.csv', "w") 

columnTitleRow = "p, s, k\n"
csv.write(columnTitleRow)

for key in b:
	p = str(key[0])
	s = str(key[1])
	k = str(key[2])
	row = p  + "," + s + "," + k + "\n"
	csv.write(row)

csv = open('beg_mid_end_syn_all1_p.csv', "w") 

columnTitleRow = "p, s, k\n"
csv.write(columnTitleRow)

for key in c:
	p = str(key[0])
	s = str(key[1])
	k = str(key[2])
	row = p  + "," + s + "," + k + "\n"
	csv.write(row)

csv = open('beg_mid_end_noalt_all1_p.csv', "w") 

columnTitleRow = "p, s, k\n"
csv.write(columnTitleRow)

for key in d:
	p = str(key[0])
	s = str(key[1])
	k = str(key[2])
	row = p  + "," + s + "," + k + "\n"
	csv.write(row)

csv = open('beg_mid_end_double_all1_p.csv', "w") 

columnTitleRow = "p, s, k\n"
csv.write(columnTitleRow)

for key in e:
	p = str(key[0])
	s = str(key[1])
	k = str(key[2])
	row = p  + "," + s + "," + k + "\n"
	csv.write(row)



csv.close()







                                                                         


                                                                         


