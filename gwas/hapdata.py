# This file describes a class for storing and manipulating haplotype data

from __future__ import division
from bisect import *
import sys
import numpy as np
from gwas import missing, complete_cases, is_na

zero=np.int16(0)
un=np.int16(1)
deux=np.int16(2)

class HapDataset():
    def __init__(self,fileName,nsnp,nhaplo):
        self.name=fileName
        ## Haplotypes
        self.nhaplo=nhaplo
        self.ihaplo=0
	## Individuals
	self.indiv={}
	self.indivIdx={} # each individual can be associated with a list of haplotypes
        ## SNPs
        self.nsnp=nsnp
        self.isnp=0
        self.snp={}
	self.snp_pos=[] # only used for fast treatment of simulated samples with one single chromosome
        self.snpIdx={} # SNP indices in Data Matrix
        self._initDataMatrix()
        ## populations
        ## populations are accessed by a name (key of the dict)
        ## the value is a vector of booleans taking a True value for haplotypes
        ## belonging to the pop
        self.populations={}
    def _initDataMatrix(self):
        ''' Initialize Data Matrix to an int16 array filled with -1'''
        self.Data=-np.ones(shape=(self.nhaplo,self.nsnp),dtype='int16')
    def addSnp(self,Name):
        ''' Add a snp to the dataset '''
        if self.snp.has_key(Name):
            return
        self.snp[Name]=SNP(Name)
        self.snpIdx[Name]=self.isnp
        self.isnp+=1
    def rmSnp(self,Name):
        ''' 
        Remove a snp from the dataset 
        The SNP is made unaccessible but the data stay in the Data Matrix
        '''
        if not self.snp.has_key(Name):
            return
        else:
            del self.snp[Name]
            del self.snpIdx[Name]
            del self._allSnps
    def SetHaplotype(self,i,snp,char_all):
        '''
        Set haplotype number i at snp SNP to value VALUE
        '''
        try:
            j=self.snpIdx[snp]
            self.Data[i,j]=self.snp[snp].addHapObs(char_all)     
        except KeyError:
            print 'snp or haplotype does not exist'
            raise
    def SetHaplotype_fast(self,i,j,snp,char_all):
        '''
        Set haplotype number i at snp j to value VALUE
	Assumes that we know the name of snp j
        '''
        try:
            self.Data[i,j]=self.snp[snp].addHapObs(char_all)     
        except KeyError:
            print 'snp or haplotype does not exist'
            raise
    def GetHaplotype(self,i,snp):
        '''
        Returns the value of haplotype i at a given SNP 
        '''
        try:
            j=self.snpIdx[snp]
	    try:
		return self.Data[i,j]<0 and missing or self.Data[i,j]
	    except KeyError:
		print 'Haplotype ',i,' not found'
		return missing
        except KeyError:
	    print 'SNP ',snp,'not found'
            return missing
    def GetHaplotypes(self,snp,pop=None):
        '''
        Return the value of all haplotypes in the dataset (or in a population) at a given SNP
        '''
	if pop is None:        
	    try:
                j=self.snpIdx[snp]
		return [ (x<0 and missing) or x for x in self.Data[:,j]]
            except KeyError:
		print 'SNP ',snp,'not found'
                return [missing]*self.nhaplo           
	else:
	    try:
		vpop=self.populations[pop]
	    except KeyError:
		print 'invalid population required'
		return
	    try:
                j=self.snpIdx[snp]
            except KeyError:
		print 'SNP ',snp,'not found'
                return [missing]*sum(vpop)
            return [ (x<0 and missing) or x for x in self.Data[vpop,j]]
    def UpdatePop(self,i,pop):
	'''
        Attach haplotype number i to a given population
        '''
	try:
            vpop=self.populations[pop]
        except KeyError:
            vpop=np.array([0]*self.nhaplo,dtype=bool)
            self.populations[pop]=vpop
        vpop[i]=True
    def addIndividual(self,ID,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        ''' Add an individual to the dataset'''
        if self.indiv.has_key(ID):
            return
        self.indiv[ID]=Individual(ID,pop,sex,fatherID,motherID,phenotype)
        self.indivIdx[ID]=[self.ihaplo,self.ihaplo+1]
        if pop is not None:
            self.UpdatePop(self.ihaplo,pop)
	    self.UpdatePop(self.ihaplo+1,pop)
	self.ihaplo+=2
    def updateIndividual(self,ID,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        ''' update individual propoerties (other than genotype)'''
	try:
	    self.indiv[ID].pop=pop
            self.indiv[ID].sex=sex
	    self.indiv[ID].fatherID=fatherID
	    self.indiv[ID].motherID=motherID
	    self.indiv[ID].phenotype=phenotype
            if pop is not None:
                haplos=self.indivIdx[ID]
                self.UpdatePop(haplos[0],pop)
	        self.UpdatePop(haplos[1],pop)
	except KeyError:
	    print 'trying to update an individual which does not exist'
    def SetIndivHaplotypes(self,ID,snp,char_all1,char_all2):
        '''
        Set the haplotypes of one individual at a snp 
        '''
        try:
            j=self.snpIdx[snp]
	    haplos=self.indivIdx[ID]
            self.Data[haplos[0],j]=self.snp[snp].addHapObs(char_all1)
	    self.Data[haplos[1],j]=self.snp[snp].addHapObs(char_all2)   
        except KeyError:
            print 'snp or individual does not exist'
            raise
    def GetIndivHaplotypes(self,ID,snp):
        '''
        Returns the two haplotypes of one individual at a snp
        '''
	try:
            j=self.snpIdx[snp]
	    try:
		haplos=self.indivIdx[IO]
		return [self.Data[haplos[0],j]<0 and missing or self.Data[haplos[0],j],self.Data[haplos[1],j]<0 and missing or self.Data[haplos[1],j]]
	    except KeyError:
		print 'Individual ',ID,' not found'
		return [missing,missing]
        except KeyError:
	    print 'SNP ',snp,'not found'
            return [missing,missing]
    def introduce_errors(self,prob_err):
	'''
	introduces random errors with probabiliy prob_err
	'''
	err=np.random.binomial(1,prob_err,size=(self.nhaplo,self.nsnp))
	self.Data=(1-self.Data)*err+self.Data*(1-err)

# population genetics routines

def Haplos_and_counts(dataset,map,pop_name,mac=1):
    '''
    For an object of class hapdata and a related map, this function returns for each chromosome :
	- a 1D array of the polymorphic positions with minor allele count greater than mac, sorted by physical position.
	- a 1D array of minor allele counts for the same positions
	- a 2D array with all haplotypes in the data set, for the same positions.
    As there are several chromosomes, outputs are lists of arrays.  
    '''
    try:
	vpop=dataset.populations[pop_name]
    except KeyError:
	print 'invalid population required'
    n=sum(vpop)
    v_counts=np.sum(dataset.Data[vpop,],axis=0)
    w=complete_cases(dataset.Data[vpop,])
    sel_snp=np.prod(w,axis=0,dtype='bool')*(v_counts>=mac)*(v_counts<=n-mac)
    u=n-v_counts
    w=(v_counts>u)
    v_counts[w]=u[w]
    # give a list of counts, haps and positions for each chromosome (sorted by position)
    pos_list=[]
    count_list=[]
    hap_list=[]
    for chro in map.chromosomes:
	snp_idx=[]
	pos=[]
	for snp in map.physMap(chro):
	    j=dataset.snpIdx[snp]
	    if sel_snp[j]:
		snp_idx.append(j)
		marker=map.position(snp)
		pos.append(marker[2])
    	pos_list.append(np.array(pos,dtype='int32'))
	count_list.append(v_counts[snp_idx])
	hap_list.append(dataset.Data[:,snp_idx][vpop,:])
    return pos_list,count_list,hap_list

def Haplos_and_counts_fast(dataset,pop_name,mac=1,prob_error=0):
    '''
    Fast version of Haplos_and_counts. It assumes that :
	- the dataset includes one single chromosome
	- SNP in the dataset are already sorted by position, and are called snp1,snp2 ...
    	- all SNP have complete information (no missing data)
    As there is a single chromosome, the outputs are arrays and not lists of arrays.
    '''
    try:
	vpop=dataset.populations[pop_name]
    except KeyError:
	print 'invalid population required'
    n=sum(vpop)
    counts=np.sum(dataset.Data[vpop,],axis=0)
    sel_snp=(counts>=mac)*(counts<=n-mac)
    u=n-counts
    w=(counts>u)
    counts[w]=u[w]
    haps=dataset.Data[vpop,:]
    return dataset.snp_pos[sel_snp],counts[sel_snp],haps[:,sel_snp]

def unfolded_counts_fast(dataset,pop_name,mac=1):
    '''
    computes unfolded counts, assuming that :
	- dataset includes one single chromosome
    	- all snps have complete information
    does not return lists of arrays but single arrays.
    '''
    try:
	vpop=dataset.populations[pop_name]
    except KeyError:
	print 'invalid population required'
    n=sum(vpop)
    counts=np.sum(dataset.Data[vpop,],axis=0)
    sel_snp=(counts>=mac)*(counts<=n-mac)
    return counts[sel_snp]

def counts_fast(dataset,pop_name,mac=1):
    '''
    computes folded counts, assuming that :
	- dataset includes one single chromosome
    	- all snps have complete information
    does not return lists of arrays but single arrays.
    '''
    try:
	vpop=dataset.populations[pop_name]
    except KeyError:
	print 'invalid population required'
    n=sum(vpop)
    counts=np.sum(dataset.Data[vpop,],axis=0)
    sel_snp=(counts>=mac)*(counts<=n-mac)
    u=n-counts
    w=(counts>u)
    counts[w]=u[w]
    return counts[sel_snp]

class SNP():
    def __init__(self,Name,chr=None,pos=None):
        self.name=Name
        self.alleles=[None,None] ## maps 0,1 to alleles
        self.mapping={} ## maps allele to 0,1
        self.addHapObs=self.addHapObsNotFull
	self.ancestral=missing ## ancestral allele - provided as a number
    def initAlleles(self,a1,a2):
        self.alleles[0]=a1
        self.alleles[1]=a2
        self.mapping[a1]=zero
        self.mapping[a2]=un
        self.addHapObs=self.addHapObsFull 
    def addHapObsFull(self,char_all):
        try:
            return self.mapping[char_all]
        except KeyError:
            return missing
    def addHapObsNotFull(self,char_all):
        ''' 
        Add an observed allele ('A','1'...) to the SNP 
        Return the 'integer haplotype',0,1
        if two alleles have already been observed, switch to fast
        version
        '''
        a=missing
        try:
            a=self.mapping[char_all]
        except KeyError:
            if self.alleles[0] is None:
                self.alleles[0]=char_all
                self.mapping[char_all]=zero
                a=zero
            elif self.alleles[1] is None:
                self.alleles[1]=char_all
                self.mapping[char_all]=un
                a=un
                self.addHapObs=self.addHapObsFull
            else:
                return a
    def addHapObsFast(self,int_all):
	'''
	Return the intger haplotype 0,1. Does not check that this value is correct (e.g. 0 or 1)
	'''
        return int16(int_all)

class Locus():
    def __init__(self,nom,position):
        self.name=nom
        self.pos=position
    def __cmp__(self,other):
        return cmp(self.pos,other.pos)

class Map():
    def __init__(self):
        '''
        Map class contains maps (duh)
        The chromosome names are stored in the chromosomes attribute (in sorted order)
        Each chromosome has two maps (genetic/RH and physical)
        '''
        self.chromosomes=[] 
        self.cartes={}
        self.markers={}
    def write(self,stream=sys.stdout):
        for chrom in self.chromosomes:
            for mk in self.physMap(chrom):
                pos=self.position(mk)
                print >>stream,'\t'.join([str(x) for x in [pos[0],mk,int(pos[1]),int(pos[2])]])
        return      
    def addMarker(self,M,C,posG=0,posP=0):
        '''
        Add a marker on chromosome C, with genetic position posG and
        physical position posP
        '''
        if not self.cartes.has_key(C):
            insort(self.chromosomes,C)
            self.cartes[C]={'Genet':[],'Phys':[]}
        locG=Locus(M,posG)
        locM=Locus(M,posP)
        self.markers[M]=(C,locG,locM)
        insort(self.cartes[C]['Genet'],locG)
        insort(self.cartes[C]['Phys'],locM)
    def delMarker(self,M):
        try:
            locus = self.markers[M]
        except KeyError:
            print "No such marker"
        del self.markers[M]
        self.cartes[locus[0]]['Genet'].remove(locus[1])
        self.cartes[locus[0]]['Phys'].remove(locus[2])          
    def genetMap(self,C,xleft=-1,xright=np.inf):
        '''
        Returns and iterator on the genetic map between position xleft and xright 
        on chromosome C
        '''
        try:
            carte=self.cartes[C]['Genet']
        except KeyError:
            return []
        mk=carte[bisect_left(carte,Locus(0,xleft)):bisect_right(carte,Locus(0,xright))]
        return (m.name for m in mk)
    def physMap(self,C,xleft=0,xright=np.inf):
        '''
        Returns and iterator on the genetic map between position xleft and xright 
        on chromosome C
        '''
        try:
            carte=self.cartes[C]['Phys']
        except KeyError:
            return []
        mk=carte[bisect_left(carte,Locus(0,xleft)):bisect_right(carte,Locus(0,xright))]
        return (m.name for m in mk)
    def length(self,C):
        '''
        Return the length of maps of chromosome C
        '''
        carte=self.cartes[C]['Phys']
        lenP=carte[-1].pos-carte[0].pos
        carte=self.cartes[C]['Genet']
        lenG=carte[-1].pos-carte[0].pos
        return [lenP,lenG]
    def chromosomes(self):
        return self.cartes.key()
    def position(self,M):
        '''
        Returns chromosome,genetic,physical position of marker M
        '''
        try:
            locus=self.markers[M]
            return locus[0],locus[1].pos,locus[2].pos
        except KeyError:
            print M,'not found'
            return None

class Individual():
    def __init__(self,id,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        self.id=id
        self.pop=pop
        self.sex=sex
        self.fatherID=fatherID
        self.motherID=motherID
        self.phenotype=[phenotype]
        self.callrate=0

