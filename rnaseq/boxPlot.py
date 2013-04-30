import pylab as pl
import os,sys
import math

def getData(filename,targetAnn):
	f = open(filename)
	data = []
	specific = []
	nonspecific = []
	foldspecific = []
	foldnonspecific = []
	lengths = []
	for i in targetAnn:
		data.append([])
		specific.append([])
		nonspecific.append([])
		foldspecific.append([])
		foldnonspecific.append([])
		lengths.append([])
	for r in f:
		if r.startswith('#'):
			continue
		tokens = r.strip().split('\t')
		if 'RRNA' in tokens[-1].upper() or 'TRNA' in tokens[-1].upper():
			continue
		else:
			ann = tokens[-2]
			if ann.isdigit():
				ann = tokens[-1]
			for i,ta in enumerate(targetAnn):
				if ann.upper() == ta.upper():
					if int(tokens[4]) >= 10:
						logv = math.log(int(tokens[4]))
						data[i].append(logv)
						lengths[i].append(int(tokens[2])-int(tokens[1]))
						if int(tokens[7]) < 10 and int(tokens[8])<10:
							specific[i].append(logv)
						else:
							nonspecific[i].append(logv)
						if (int(tokens[7]) == 0 or int(tokens[4]) *1.0 / int(tokens[7]) >= 2) and (int(tokens[8]) ==0 or int(tokens[4])*1.0/int(tokens[8]) >= 2):
							foldspecific[i].append(logv)
						else:
							foldnonspecific[i].append(logv)
					break
	f.close()
	return data,specific,nonspecific,foldspecific,foldnonspecific,lengths
			
			
def plotSpec(specific,nonspecific,prefixes,categories,stype):
	for index,c in enumerate(categories):
		fig = pl.figure()
		if stype == 'abs':
			pl.title((c + ' threshold').title())
		else:
			pl.title((c+' fold').title())
		pl.boxplot(specific[index]+nonspecific[index])
		print c, ' ', stype
		for j,pre in enumerate(prefixes):
			print '\t', pre,'\t',len(specific[index][j]),len(nonspecific[index][j])
		pl.ylim([2,15])
		pl.ylabel('Log Expression Level')
		pl.xticks([i+1 for i in range(len(specific[index]+nonspecific[index]))],[i+' Specific' for i in prefixes] + [i+' Non-specific' for i in prefixes],rotation=30)
		fig.autofmt_xdate()
		fig.savefig(c+'.'+stype+'.spec.png',dpi=fig.dpi)
	

def main():
	#exon, intron, unknown
	specific = [[],[],[]]
	nonspecific = [[],[],[]]
	foldspecific = [[],[],[]]
	foldnonspecific = [[],[],[]]
	for prefix in sys.argv[1:]:
		tempexonic,tempspecific,tempnonspecific,tempfoldspecific,tempfoldnonspecific,templength = getData(prefix+'.exonic.overlap.out.annotation.txt',['exon',])
		exonic = tempexonic[0]
		exoniclength = templength[0]
		specific[0].append(tempspecific[0])
		nonspecific[0].append(tempnonspecific[0])
		foldspecific[0].append(tempfoldspecific[0])
		foldnonspecific[0].append(tempfoldnonspecific[0])
		tempData,tempspecific,tempnonspecific,tempfoldspecific,tempfoldnonspecific,templength = getData(prefix + '.novel.overlap.out.annotation.txt',['intron','unknown'])
		intronic = tempData[0]
		unknown = tempData[1]
		introniclength = templength[0]
		unknownlength = templength[1]
		specific[1].append(tempspecific[0])
		specific[2].append(tempspecific[1])
		nonspecific[1].append(tempnonspecific[0])
		nonspecific[2].append(tempnonspecific[1])
		foldspecific[1].append(tempfoldspecific[0])
		foldspecific[2].append(tempfoldspecific[1])
		foldnonspecific[1].append(tempfoldnonspecific[0])
		foldnonspecific[2].append(tempfoldnonspecific[1])
		plotData = [exonic,intronic,unknown]
		print prefix
		print "exonic: ", len(exonic)
		print 'intronic: ', len(intronic)
		print 'unknown: ', len(unknown)
		fig = pl.figure()
		pl.boxplot(plotData)
		pl.ylim([2,15])
		pl.ylabel('Log Expression Level')
		pl.xticks([1,2,3],['Exonic','Intronic','Unknown'])
		pl.title(prefix.replace('_fsorted','').replace('_',' '))
		fig.savefig(prefix+'.expression.png',dpi=fig.dpi)
		
		fig = pl.figure()
		pl.boxplot([exoniclength,introniclength,unknownlength])
		pl.ylim([60,2500])
		pl.ylabel('Transcript Length')
		pl.xticks([1,2,3],['Exonic','Intronic','Unknown'])
		pl.title(prefix.replace('_fsorted','').replace('_',' '))
		fig.savefig(prefix+'.length.png',dpi=fig.dpi)
		#pl.show()
	abbr = []
	for i in sys.argv[1:]:
		tokens = i.split('_')
		abbr.append(tokens[0][0].upper() + tokens[1][0:2].title())
	plotSpec(specific,nonspecific,abbr, ['exonic','intronic','unknown'],'abs')
	plotSpec(foldspecific,foldnonspecific,abbr, ['exonic','intronic','unknown'],'fold')
			
	

if __name__=='__main__':
	main()

