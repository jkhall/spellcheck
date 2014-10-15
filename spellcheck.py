import sys
import re
import argparse
import string
import math

"""
I'm controlling the number of lines I'm reading from the file, within
computeCorrectSpelling, and I have it set to 100

Tuning:

I tried to account for close keys in the matrix such as, a-s, g-h, and s-d

I lowered the costs on vowel to vowel substitutions

I also added a weight parameter, but i could not determine a way to get it 
to improve my results

hypothesis.txt is the first one hundred results for my basic bigram

hypothesis2.txt is the first hundred for bigram with my altered confusion matrix

"""


#metainfo
__author__= "Jordan Hall"
__date__= "Oct 2014"
__version__= "1"

def buildSModel(n, sourcefile):
	#read in doc
	#build out dictionary with probs doing either unigram or bigran
	foundstart = False;
	hit2gram = False;	
	sourcedict = {}
	wordcount = 0	
	if n == 1:
		count = 0;
		for line in open(sourcefile):
			if line.strip("\\") == "2-grams:\n":
				print line
				break
			if foundstart:	
				linelist = line.split();
				if len(linelist) == 1 or linelist == []:
					print "odd line"
					continue
				sourcedict[linelist[1]] = float(linelist[0]);
			else:
			#matchObjOne = re.match('\\2-grams:', line);
			#matchObjTwo = re.match('\\1-grams:', line);
				if not line.strip("\\") == "1-grams:\n": #think that might be an escape char
					continue #go to next iteration
				else: #you have the correct start line, build
					foundstart = True;
		return sourcedict;	
	if n == 2:
		#have to keep track of bow, and conditional prob
		for line in open(sourcefile):
			if line.strip("\\") == "1-grams:\n":
				foundstart = not foundstart;
				continue
			if line.strip("\\") == "2-grams:\n":
				foundstart = not foundstart;
				hit2gram = not hit2gram;
				print "hit 2gram"
				continue;
			if foundstart:
				linelist = line.split();
				if len(linelist) == 1 or linelist == []:
					continue
				try:
					sourcedict[linelist[1]]
				except KeyError:
					sourcedict[linelist[1]] = {}
				if len(linelist) < 3:
					sourcedict[linelist[1]]["**n_prob_val**"] = float(linelist[0]); #still need this in case of backoff
					continue;
				sourcedict[linelist[1]]["**bow_val**"] = float(linelist[2]);
				sourcedict[linelist[1]]["**n_prob_val**"] = float(linelist[0]); #still need this in case of backoff
			if hit2gram:
				#arpa format is prevword word
				#sourcedict[prevword][word]
				linelist = line.split()
				if len(linelist) == 1 or linelist == []:
					continue
				try:
					sourcedict[linelist[1]]
				except KeyError:
					sourcedict[linelist[1]] = {}
				sourcedict[linelist[1]][linelist[2]] = float(linelist[0]);

		return sourcedict;

def buildCModel(channelfile):
	#iterate over the two words and compute which one is the correct one
	edistdict = {letter:{} for letter in string.lowercase}
	edistdict["eps"] = {}
	lineskip = True
	for line in open(channelfile):
		if lineskip:
			elementlist = line.split();
			lineskip = False;
			continue
		linelist = line.split();
		for i in range(len(elementlist)):
			edistdict[linelist[0]][elementlist[i]] = float(linelist[i + 1]); #just to skip that leading 0
	return edistdict;		
				
def computeEditDist(errorword, correctword, cmodel):
	#iterate through model computing edit distance for the word
	if not errorword.isalpha():
		return
	if not correctword.isalpha():
		return
	cwordlist = ["eps"]
	cwordlist += list(correctword);
	ewordlist = ["eps"]
	ewordlist+= list(errorword);
	editdist = 0.00;
	#initvals = range(len(ewordlist))
	#toprowvals = range(len(cwordlist))
	#currowvals = range(len(cwordlist))
	mtrix = range(len(ewordlist))
	#would like to avoid all of this looping if at all possible
	for row in range(len(ewordlist)):
		mtrix[row] = range(len(cwordlist));
	mtrix[0][0] = 0; #to avoid log issue
	#should probably have various edit distance probs
	#for edge cases	
	#if len(ewordlist) == min(len(ewordlist), len(cwordlist)):
		#fill with nulls
		#while not len(ewordlist) == len(cwordlist):
			#ewordlist.append("eps");
	#else:
		#fill the other list with nulls
		#while not len(ewordlist) == len(cwordlist):
			#cwordlist.append("eps")
	#append null to start
		
	for i in range(1, len(ewordlist)):
		#compute P(w | c)
		#need to actually compute edit dist
		#build init vals
		#don't forget to add backwards values
		mtrix[i][0] = math.log(cmodel[cwordlist[0]][ewordlist[i]]) + mtrix[i - 1][0];
		#initvals[i] = -math.log(float(cmodel[ewordlist[i]][cwordlist[0]])) + initvals[i-1]; #+= bc they add up
		#else:
			#initvals[i] = -math.log(float(cmodel[ewordlist[i]][cwordlist[0]])); #+= bc they add up
	for i in range(1, len(cwordlist)):	
		#if i == 0:
		#	mtrix[0][i] = math.log(cmodel[cwordlist[i]][ewordlist[0]])
		#else: 
		mtrix[0][i] = math.log(cmodel[cwordlist[i]][ewordlist[0]]) + mtrix[0][i-1]
		#if not i == 0:
		#toprowvals[i] = -math.log(float(cmodel[ewordlist[0]][cwordlist[i]])) + toprowvals[i - 1]; 
		#else:
		#	toprowvals[i] = -math.log(float(cmodel[ewordlist[0]][cwordlist[i]])); 
	#you want to travel down init vals, computing at each step, but across the c word first
	for j in range(1, len(ewordlist)):
		#currowvals[0] = initvals[j];
		# have to consider currowvals[i-1], initvals[ewordindex-1], toprowvals[i]
		#at some point you have to know the costs for each of the ops
		for k in range(1, len(cwordlist)):
			#should you be able to compute the costs once? I think so.
			costlist = []
			delcost = math.log(cmodel[cwordlist[k]]["eps"])
			costlist.append(delcost)
			inscost = math.log(cmodel["eps"]["k"]) #needing some letter
			costlist.append(inscost)
			submatchcost = math.log(cmodel[ewordlist[j]][cwordlist[k]]);
			#except KeyError:
			#	if ewordlist[j] == cwordlist[k]:
			#		submatchcost = -math.log(0.75);
			#	else:
			#		submatchcost = -math.log(0.005);
			costlist.append(submatchcost)
			costlist[0] += mtrix[j][k-1] 
			costlist[1] += mtrix[j-1][k]
			costlist[2] += mtrix[j-1][k-1]
			mtrix[j][k] = max(costlist); #pick the minimum cost
			#print "Pick:" + str(mtrix[j][k])	
			#do the above until you have that final cell and return that as the editdist
			#switch currowvals and toprow
		#toprowvals = currowvals;
	
		#try:
		#	editdist += math.log(float(cmodel[ewordlist[i]][cwordlist[i]]));
		#except KeyError:
			#add deletion cost
	editdist = mtrix[-1][-1]; #last item in the list?
	#print "Edit dist: " + str(editdist)
	return editdist

def computeCorrectSpelling(correctfile, outfile, cmodel, smodel, n, weight=1):
	currentcalc = 0;
	correctspelling = ""
	errorstring = '<ERROR>'
	f = open(outfile, 'w+')
	count = 0
	errordict = {}	
	for line in open("penntreebank.errors/" + correctfile):
		if count == 100: break;
		count += 1
		cline = line.split();
		#with line split find '<Error>'
		for i in range(len(cline)):
			if set(errorstring).issubset(set(cline[i])):
				#extract word from within
				#errorword = cline[i+1]
				#build correct string
				errorword = re.findall(">.*<", cline[i])
				errorword = errorword[0].strip("><")
				high = -1000
				if n == 1:
					try:
						errordict[errorword]
					except KeyError:
						for correctword in smodel.keys():
							if correctword == errorword: #base assumption that the word is wrong
								continue
							currentED = computeEditDist(errorword, correctword, cmodel)
							currentPC = smodel[correctword]
							if currentED is None:
								continue
							currentcalc = weight*currentED + currentPC; 	
							if currentcalc > high: 
								high = currentcalc; 
								correctspelling = correctword;
						errordict[errorword] = correctspelling;
					correctspelling = errordict[errorword]
					f.write(correctspelling + "\n")
				if n == 2:
					#use prevword and nextword in consideration
					high = -1000
					try:
						errordict[errorword]
					except KeyError:
						for potentialword in smodel.keys():
							if potentialword == errorword:
								continue;
							try:
								prevword = cline[i-1]#assuming there is a prev word
							except IndexError:
								prevword = "<unk>"
							try:
								nextword = cline[i+1]#what if this is an error?
							except IndexError:
								nextword = "<unk>"
							currentED = computeEditDist(errorword, potentialword, cmodel);
							if currentED is None:
								continue
							currentPC = None
							currentPrevP = None
							currentNextP = None
							#try hash using prevword -> if fail use backoff and unk
							try:
								smodel[prevword][potentialword]
							except KeyError:
								#word
								try: 
									smodel[prevword]
								except KeyError:
									currentPrevP = smodel["<unk>"]["**n_prob_val**"]+smodel[potentialword]["**n_prob_val**"]
								if currentPrevP is None:
									currentPrevP = smodel[prevword]["**bow_val**"]+smodel[potentialword]["**n_prob_val**"]
							#try hash using nextword -> if fail use backoff and unk
							try: 
								smodel[potentialword][nextword] #<-is it the context or the word
							except KeyError:
								#word
								try:
									smodel[nextword];
								except KeyError:
									currentNextP = smodel["<unk>"]["**n_prob_val**"]+smodel[potentialword]["**bow_val**"]
								if currentNextP is None:
									#context prob
									currentNextP = smodel[potentialword]["**bow_val**"]+smodel[nextword]["**n_prob_val**"]	
							#try hash using sdict[prevword][potentialword] -> if fail use backoff
							#try hash using sdict[potentialword][nextword] -> if fail use backoff
							if currentNextP is None:
								#print "yes"
								currentNextP = smodel[potentialword][nextword]
							if currentPrevP is None:
								#print "yes"
								currentPrevP = smodel[prevword][potentialword]
							currentPC = currentPrevP + currentNextP
							currentcalc = currentPC + currentED;
							if currentcalc > high:
								high = currentcalc
								correctspelling = potentialword;		 
						errordict[errorword] = correctspelling
					correctspelling = errordict[errorword];
					f.write(correctspelling + "\n")
#really all the channel model needs to do is compute the costs
def editdistcomp2(wrongword, rightword, channelmod):
	wronglist = ["eps"]
	wronglist += list(wrongword)
	rightlist = ["eps"]
	rightlist += list(rightword)
	matrix = []
	matrix.append(0);
		#deletion cost
			
	
if __name__ == '__main__':
	myparser = argparse.ArgumentParser(description='Command line arg parser');
	myparser.add_argument('-lmfile', required=True);
	myparser.add_argument('-n', required=True);
	myparser.add_argument('-channel', required=True);
	myparser.add_argument('-infile', required=True);
	myparser.add_argument('-o', required=True);
	myparser.add_argument('-w'); #optional
	args = vars(myparser.parse_args());
	sourcefile = args['lmfile'];
	channelfile = args['channel'];
	correctfile = args['infile']
	outfile = args['o'];
	n = int(args['n']);
	try:
		weight = float(args['w'])
	except TypeError:
		weight = None #think this ought to be suitable
	sourcedict = buildSModel(n, sourcefile);
	channeldict = buildCModel(channelfile);
	#print channeldict["eps"]["k"]
	#computeEditDist("th", "the", channeldict)
	#print sourcedict["to"]
	#print sourcedict["government"]
	if not weight is None:
		computeCorrectSpelling(correctfile, outfile, channeldict, sourcedict, n, weight);
	else:
		computeCorrectSpelling(correctfile, outfile, channeldict, sourcedict, n); #default val is 1
