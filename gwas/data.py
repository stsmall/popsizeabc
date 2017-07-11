# This file describes a class for storing and manipulating genotype data

from __future__ import division
from bisect import *
import sys
import numpy as np
from gwas import missing, complete_cases, is_na

zero=np.int16(0)
un=np.int16(1)
deux=np.int16(2)

class Dataset():
    def __init__(self,fileName,nsnp,nindiv):
        self.name=fileName
        ## Individuals
        self.nindiv=nindiv
        self.iindiv=0
        self.indiv={}
	## sorted list of individuals
	self.indiv_list=[]
        ## Individual index in Data Matrix
        self.indivIdx={}
        ## SNPs
        self.nsnp=nsnp
        self.isnp=0
        self.snp={}
        ## SNP indices in Data Matrix
        self.snpIdx={}
        self._initDataMatrix()
        ## populations
        ## populations are accessed by a name (key of the dict)
        ## the value is a vector of booleans taking a True value for individuals
        ## belonging to the pop otherwise
        ## These vectors can be used for fast access to population genotypes e.g.:
        ## vpop=self.populations['mypop']
        ## mypop_genotypes=self.Data[vpop,]
        self.populations={}
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
    def addIndividual(self,ID,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        ''' Add an individual to the dataset'''
        if self.indiv.has_key(ID):
            return
        self.indiv[ID]=Individual(ID,pop,sex,fatherID,motherID,phenotype)
        self.indivIdx[ID]=self.iindiv
	self.indiv_list.append(ID)
        if pop is not None:
            try:
                vpop=self.populations[pop]
            except KeyError:
                vpop=np.array([0]*self.nindiv,dtype=bool)
                self.populations[pop]=vpop
            vpop[self.iindiv]=True
	self.iindiv+=1
    def updateIndividual(self,ID,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        ''' update individual propoerties (other than genotype)'''
	try:
	    self.indiv[ID].pop=pop
            self.indiv[ID].sex=sex
	    self.indiv[ID].fatherID=fatherID
	    self.indiv[ID].motherID=motherID
	    self.indiv[ID].phenotype=phenotype
            if pop is not None:
                try:
                    vpop=self.populations[pop]
                except KeyError:
                    vpop=np.array([0]*self.nindiv,dtype=bool)
                    self.populations[pop]=vpop
                vpop[self.indivIdx[ID]]=True
	except KeyError:
	    print 'trying to update an individual which does not exist'
    def _initDataMatrix(self):
        ''' Initialize Data Matrix to an int16 array filled with -1'''
        self.Data=-np.ones(shape=(self.nindiv,self.nsnp),dtype='int16')
    def IndivGenotype(self,indiv,snps=None):
        '''
        Return (0,1,2) genotype of individual INDIV at all SNPs in the dataset
        '''
        try:
            i=self.indivIdx[indiv]
        except KeyError:
            print 'Individual',indiv,'not in Data'
            raise
        if snps is None:
            return [ (x<0 and missing) or x for x in  self.Data[i,]]
        else:
            return [ (x<0 and missing) or x for x in self.Data[i,[self.snpIdx[s] for s in snps]]]
    def IndivGenotype_char(self,indiv,snps=None):
        igeno=self.IndivGenotype(indiv,snps)
        cgeno=[]
        if snps is None:
            snps=self.AllSnps()
        for i,s in enumerate(snps):
            cgeno.append(self.Alleles(igeno[i],s))
        return cgeno
    def AllSnps(self):
        try:
            return self._allSnps
        except AttributeError:
            self._allSnps=sorted(self.snp,key=lambda x: self.snpIdx[x])
            return self._allSnps
    def SnpGenotype(self,snp,pop=None):
        '''
        Return genotype of snp SNP for all individuals in the dataset or in a population
        '''
	if pop is None:        
	    try:
                j=self.snpIdx[snp]
            except KeyError:
                return [missing]*self.nindiv
            return [ (x<0 and missing) or x for x in self.Data[:,j]]
	else:
	    try:
		vpop=self.populations[pop]
	    except KeyError:
		print 'invalid population required'
		return
	    try:
                j=self.snpIdx[snp]
            except KeyError:
                return [missing]*sum(vpop)
            return [ (x<0 and missing) or x for x in self.Data[vpop,j]]
    def Alleles(self,geno,snp):
        '''
        Return the alleles corresponding to genotype GENO at snp SNP
        '''
        if geno is missing:
            return ('?','?')
        snpAll=self.snp[snp].alleles
        if geno==0:
            return (snpAll[0],snpAll[0])
        elif geno==1:
            return (snpAll[0],snpAll[1])
        else:
            assert geno==2
            return (snpAll[1],snpAll[1])
    def Genotype(self,indiv,snp):
        '''
        Return the genotype of individual INDIV  at snp SNP 
        If SNP is not found in the dataset, returns numpy.ma.masked
        '''
        try:
            i=self.indivIdx[indiv]
        except KeyError:
            print 'Individual',indiv,'not in Data'
            raise
        try:
            j=self.snpIdx[snp]
        except KeyError:
            return missing
        return self.Data[i,j]<0 and missing or self.Data[i,j]
    def SetGenotypeSlow(self,indiv,snp,value):
        '''
        Set the genotype of individual INDIV at snp SNP to value VALUE
        Correct Genotype is a value of length 2, with alleles matching the SNP
        '''
        try:
            i=self.indivIdx[indiv]
        except KeyError:
            print 'Individual',indiv,'not in Data'
            raise
        try:
            j=self.snpIdx[snp]
        except KeyError:
            print 'SNP',snp,'not in Data'
            raise
        if self.Data[i,j] > -1:
            self.snp[snp].delGenotype(self.Data[i,j])
        geno=self.snp[snp].addObservation(value)
        if geno is not missing:
            self.indiv[indiv].callrate+=1
        self.Data[i,j]=geno
    def SetGenotype(self,indiv,snp,value):
        '''
        Fast Version of SetGenotype with no checking for errors 
        Use with caution
        --
        Set the genotype of individual INDIV at snp SNP to value VALUE
        Correct Genotype is a value of length 2, with alleles matching the SNP
        '''
        try:
            i=self.indivIdx[indiv]
            j=self.snpIdx[snp]
            self.Data[i,j]=self.snp[snp].addObservation(value)     
        except:
            print 'Error'
            raise
    def printSNP(self,stream=sys.stdout):
        '''
        Print SNP information
        '''
        print >>stream,'name\tA1\tA2\tcallRate'
        for sName in sorted(self.snp,key=lambda x:self.snpIdx[x]):
            geno=self.SnpGenotype(sName)
            s=self.snp[sName]
            print >>stream,s.name,s.alleles[0],s.alleles[1],np.mean([g is not missing for g in geno])
    
    ### POPULATION GENETICS methods
    def compute_pop_frq(self):
        '''
        Compute allele frequencies for each population in the dataset, considering only
        non founders.
        '''
        try:
            fvec=self.founders
        except AttributeError:
            self.BuildPedigree()
            fvec=self.founders
        pop_names=[]
        frqs=[]
        for pname,pvec in self.populations:
            pop_names.append(pname)
            pop_founders=pvec&fvec
            w=complete_cases(self.Data[pop_founders,])
            frqs.append(0.5*np.average(self.Data[pop_founders,],axis=0,weights=w))
        return {"pops":pop_names,'freqs':np.vstack(frqs)}
   
    def AltCount_Onesnp(self,pop_name,snp_name):
        '''
        returns the number of non missing alleles and the number of alternative alleles in one pop at a given snp
        the alternative allele is the one mapping to 1 in our dataset
        '''
        genotypes=self.SnpGenotype(snp_name,pop_name)
        tc=len(genotypes)*2
        ac=0
        for g in genotypes:
	    if g == missing:
	        tc-=2
	    else:
	        ac+=g
        return [tc,ac]

def All_Counts(dataset,map,pop_name):
    '''
    returns a list of minor counts and positions for each chromosome (sorted by position)
    '''
    try:
	vpop=dataset.populations[pop_name]
    except KeyError:
	print 'invalid population required'
    n=sum(vpop)*2
    v_counts=np.sum(dataset.Data[vpop,],axis=0)
    w=complete_cases(dataset.Data[vpop,])
    sel_snp=np.prod(w,axis=0,dtype='bool')*(v_counts>0)*(v_counts<n)
    nb_snp=sum(sel_snp)
    u=n-v_counts
    w=(v_counts>u)
    v_counts[w]=u[w]
    snp_names=sorted(dataset.snpIdx,key=lambda x: dataset.snpIdx[x])
    # create a dictionnary (name,count) for selected snps
    dico_counts={}
    for i in range(0,len(v_counts)):
	if sel_snp[i]:
	    dico_counts[snp_names[i]]=v_counts[i]
    # give a list of counts and positions for each chromosome (sorted by position)
    pos_list=[]
    count_list=[]
    for chro in map.chromosomes:
	pos=[]
	counts=[]
	for snp in map.physMap(chro):
	    try:
		counts.append(dico_counts[snp])
		marker=map.position(snp)
		pos.append(marker[2])
	    except KeyError:
		pass
    	pos_list.append(np.array(pos,dtype='int32'))
	count_list.append(np.array(counts,dtype='int16'))
    return pos_list,count_list

def Genos_and_counts(dataset,map,pop_name,mac=1):
    '''
    returns a list of genotypes, minor counts and positions for each chromosome (sorted by position)
    '''
    try:
	vpop=dataset.populations[pop_name]
    except KeyError:
	print 'invalid population required'
    n=sum(vpop)*2
    v_counts=np.sum(dataset.Data[vpop,],axis=0)
    w=complete_cases(dataset.Data[vpop,])
    sel_snp=np.prod(w,axis=0,dtype='bool')*(v_counts>=mac)*(v_counts<=n-mac)
    u=n-v_counts
    w=(v_counts>u)
    v_counts[w]=u[w]
    # give a list of counts, genos and positions for each chromosome (sorted by position)
    pos_list=[]
    count_list=[]
    geno_list=[]
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
	geno_list.append(dataset.Data[:,snp_idx][vpop,:])
    return pos_list,count_list,geno_list

class Individual():
    def __init__(self,id,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        self.id=id
        self.pop=pop
        self.sex=sex
        self.fatherID=fatherID
        self.motherID=motherID
        self.phenotype=[phenotype]
        self.callrate=0

class SNP():
    def __init__(self,Name):
        self.name=Name
        self.alleles=[None,None] ## maps 0,1 to alleles
        self.mapping={} ## maps allele to 0,1
        self.counts=np.zeros(3,dtype='int16') ## number(a1),number(a2),number(missing)
        self.Gcounts=[0,0,0]
        self.addObservation=self.addObservationNotFull
	self.ancestral=missing ## ancestral allele - provided as a number
    def callrate(self):
        return 0.5*(self.counts[0]+self.counts[1])/(0.5*(self.counts[0]+self.counts[1])+self.counts[2])
    def initAlleles(self,a1,a2):
        self.alleles[0]=a1
        self.alleles[1]=a2
        self.mapping[a1]=0
        self.mapping[a2]=1
        self._create_mapping()
        self.addObservation=self.addObservationFull
    def delGenotype(self,geno):
        '''
        Remove a genotype from observed alleles
        '''
        if geno == missing:
            self.counts[2] -= 1
        else:
            self.counts[0] -= 2-geno
            self.counts[1] -= geno
    def addGenotype(self,geno):
        '''
        Add a genotype to observed alleles
        '''
        if geno==missing:
            self.counts[2]-=1
        else:
            self.counts[0] += 2-geno
            self.counts[1] += geno
    def _create_mapping(self):
        self.g2int={}
        self.g2int[self.alleles[0]+self.alleles[0]]=zero
        self.g2int[self.alleles[0]+self.alleles[1]]=un
        self.g2int[self.alleles[1]+self.alleles[0]]=un
        self.g2int[self.alleles[1]+self.alleles[1]]=deux
    def addObservationFull(self,value):
        try:
            return self.g2int[value]
        except KeyError:
            return missing
    def addObservationNotFull(self,value):
        ''' 
        Add an observed genotype ('AA','01'...) to the SNP 
        Return the 'integer genotype',0,1,2
        if two alleles have already been observed, switch to fast
        version
        '''
        geno=missing
        try:
            i1,i2=value[0],value[1]
        except :
            return geno
        try:
            g1=self.mapping[i1]
        except KeyError:
            if self.alleles[0] is None:
                self.alleles[0]=i1
                self.mapping[i1]=zero
                g1=zero
            elif self.alleles[1] is None:
                self.alleles[1]=i1
                self.mapping[i1]=un
                g1=un
                self._create_mapping()
                self.addObservation=self.addObservationFull
            else:
                return geno
        try:
            g2=self.mapping[i2]
        except KeyError:
            assert self.alleles[0]!=None
            if self.alleles[1] is None:
                self.alleles[1]=i2
                self.mapping[i2]=un
                g2=un
                self._create_mapping()
                self.addObservation=self.addObservationFull
            else:
                return geno
        geno=g1+g2
        return geno

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
        self.input_order=[]
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
        self.input_order.append(M)
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
    def sort_loci(self,loci,physical=True):
        '''
        Sort loci according to their position on the map

        Parameters:
        ---------------

        loci : a list of loci names
        physical : if True sort according to physical position, o.w. use genetic  / RH position

        Return Value:
        -----------------

        list of the same loci, sorted
        '''
        new_list=[]
        for locus in loci:
            try:
                tmp=self.markers[locus]
                if physical:
                    myloc=tmp[2]
                else:
                    myloc=tmp[1]
                try:
                    mychr=int(tmp[0])
                except:
                    mychr=tmp[0]
            except KeyError:
                continue
            new_list.append((locus,mychr,myloc.pos))
        new_list=sorted(new_list,key=itemgetter(1,2))
        return [x[0] for x in new_list]
                    

