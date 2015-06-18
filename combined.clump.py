#! /usr/bin/env python
'''
cluster.py
Uses pamk from clump.py to determine clustering for mutation positions in a given gene
'''

from optparse import OptionParser
#from clump import pamk
import rpy2.robjects as robjects
import math
import sys, os
#import random
import numpy as np
r=robjects.r
r.library('fpc')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
class CmdOpts(object):
    usage="""%prog [options] -f -p


    """

    def __init__(self):
        parser = OptionParser(usage=CmdOpts.usage)
        parser.add_option("-f", "--fname", dest="filename",
                          help="""Input annotation file""")
        parser.add_option("-a","--aname",dest="affreq",
                          help="""Input allele frequency cutoff""",default=1)
        parser.add_option("-p","--pname",dest='proteinlengths',
                          help="""Input file of NP and protein lengths""")
        parser.add_option('-c','--cname',dest='controlperms',
                          help='''Input the controls that will be used for permutation testing''')
        parser.add_option('-z','--zname',dest='permutations',
                          help='''Input the number of permutations you want to perform on each gene''')
        parser.add_option("-m","--mname",dest='nmutcutoff', default=5,
                          help="""Number of mutations needed in a gene to be considered""")
        parser.add_option('-n',action="store_true",dest='normalize',
                          help="""Normalize protein lengths""")
        parser.add_option('-t',action='store_true',dest='titles',help="""Output Column Titles""")
        (opts, args) = parser.parse_args()

        if not opts.filename and not opts.affreq and not opts.proteinlengths:
            parser.error("Incorrect input args")

        self.__dict__.update(opts.__dict__)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
class Parser():
    '''
    Parses the AD_final and AR_final files
    '''
    def FillLen(self,Proteinfile):
        lendb=dict()
        for line in file(Proteinfile):
            lineinfo=line.strip().split('\t')
            lendb[lineinfo[0]]=int(lineinfo[1])
        return lendb

    def fill(self, testfile,afcutoff):
        db=dict()
        domain=dict()
        for line in file(testfile):
            
            lineinfo = line.strip().split('\t')
            key=(lineinfo[0],lineinfo[1],lineinfo[2])
            
            if len(lineinfo)==8:
                allelefreq=0
            elif lineinfo[8]=="NA":
                allelefreq=0
            else:
                allelefreq=float(lineinfo[8])
            if allelefreq<=afcutoff:
                pos = int(lineinfo[3])
                protID = lineinfo[1]
                if self.normalize:
                    plen = self.Proteinlen[protID]
                    pos=float(pos)/float(plen)
                if key in db.keys():
                    db[key].append(pos)
                    if len(lineinfo)==9:
                        if lineinfo[9]!="NA" and lineinfo[9]!="NONE":
                            domain[protID]=domain[protID]+1
                else:
                    db[key]=[pos]
                    if len(lineinfo)==9:
                        if lineinfo[9]!="NA" and lineinfo[9]!="NONE":
                            domain[protID]=1
                    else:
                        domain[protID]=0

        return db,domain

    def __init__(self, inputfile,lengthfile,af,ncutoff,normalize):
        

        self.Proteinlen = self.FillLen(lengthfile)
        self.normalize=normalize
        self.NMUTCUTOFF=ncutoff
        #Initializing dictionary 
        self.db = dict()
        #domain db
        self.domain=dict()
        #Initializing medoids
        self.medoids = dict()
        #Initializing medoids dict for easy calling (keyed by NPid)
        self.NP2medoids = dict()
        #Initializing clump scores
        self.clumps = dict()
        #Now reading in file:
        self.db,self.domain=self.fill(inputfile,af)
        #Now that the dictionary is filled, applying clustering to the mutation positions
        self.medoids = dict()
        #Grabbing the medoids for this particular key
        for key in self.db:
            sys.stderr.write('Now clustering: ' + str(key[0]) + '\n')
            muts = self.db[key]
            
            if len(muts) < int(self.NMUTCUTOFF): continue
            medoids = self.pamk(self.db[key])
            self.medoids[key] = medoids
            protID = key[1]
            self.NP2medoids[protID] = medoids

        #Applying clump
        self.clump()


    def __call__(self):
        '''
        Prints out the clump score for each record
        '''
        for record in self.clumps:
            print '\t'.join(record) + '\t' + str(self.clumps[record])

    def getmedoids(self, protID):
        try:
            return self.NP2medoids[protID]
        except KeyError:
            return None
        

    def pamk(self,data):
        '''
        Wrapper around the fpc library's pamk function '''
        data = [float(x) for x in data]
        if len(data) <= 2:
            mean = sum(data)/len(data)
            return [mean]

        r = robjects.r
        m = robjects.FloatVector(data)
        if len(data) < 4:
            krange = [2]
        else:
            kmax = int(round(len(data)/2.0))
            krange = range(2,min(kmax+1, 15))
        try:
            if len(data) > 1000:
                reply = r.pamk(m, robjects.IntVector(krange), usepam = False)
                pamobject = reply[0]
                medoids = list(pamobject[1])
            else:
                reply = r.pamk(m, robjects.IntVector(krange))
                pamobject = reply
                pamobject = pamobject[0]
                medoids = list(pamobject[0])
            return medoids
        except:
            raise
            return None


    def clump(self):
        '''
        Applies the CLUMP score for all records
        This CLUMP score is the average of the log distance of a mutation to its nearest medoid
        '''
        sys.stderr.write('Now clumping\n')
        
        for record in self.medoids:
            
            distances = []
            medoids = self.medoids[record]
            muts = self.db[record]
            for mutpos in muts:
                dist = math.log( min( [abs(mutpos - x) for x in medoids] ) + 1)
                distances.append(dist)
            #Now taking the average
            clump = float(sum(distances)) / float(len(muts))
            self.clumps[record] = clump


    def permparse(self,filehandle,afcutoff):
        '''parse file to get db for permutations'''
        db=dict()
        for line in file(filehandle):
            lineinfo = line.strip().split('\t')
            key=(lineinfo[0],lineinfo[1],lineinfo[2])
            allelefreq=float(lineinfo[8])
            if allelefreq<=afcutoff:
                protID = lineinfo[1]
                if protID in db.keys():
                    db[protID]=db[protID]+1
                else:
                    db[protID]=1
        return db

    def controlclump(self,db):
        '''perform clump on the control population'''

        clumps=dict()
        for key in db:
            medoids=self.pamk(db[key])
            distances = []
            muts = db[key]
            for mutpos in muts:
                dist = math.log( min( [abs(mutpos - x) for x in medoids] ) + 1)
                distances.append(dist)
            #Now taking the average                                                                                             
            clump = float(sum(distances)) / float(len(muts))
            clumps[key] = clump
        return clumps


    def permutedb(self,db,numberofperm):
        '''do permutations'''
        clumppermutations=dict()
        for key in db:
            if db[key]>=int(self.NMUTCUTOFF):##this was added to correct when controls have less than ncutoffs
                clumppermutations[key]=[]
                for ii in xrange(0,numberofperm):
                    try:
                    #permsample=random.sample(range(self.Proteinlen[key]),controldb[key])                                       
                        permsample=np.random.choice(range(self.Proteinlen[key]),size=db[key],replace=True)
                        permmedoids=self.pamk(permsample)
                        distances=[]
                        for perm in permsample:
                            dist = math.log( min( [abs(perm - x) for x in permmedoids] ) + 1)
                            distances.append(dist)
                        clump = float(sum(distances)) / float(len(permsample))
                        clumppermutations[key].append(clump)
                    except KeyError:
                        pass
        return clumppermutations

    def performpermutations(self,casefile,controlfile,numberofperm,afcutoff):
        '''Perform permutation testing to get a pvalue'''

        controldb,controldomains=self.fill(controlfile,afcutoff)
        controlname=controldb.keys()[0][2]
        controlvalues=self.controlclump(controldb)
        controldb=self.permparse(controlfile,afcutoff)
        casedb=self.permparse(casefile,afcutoff)
        controlpermutations=self.permutedb(controldb,numberofperm)
        casepermutations=self.permutedb(casedb,numberofperm)
        
        for record in self.clumps:
            
            protID=record[1]
            try:####added when controlpermutations[protID] gives keyerror
                controls=np.array(controlpermutations[protID])
                cases=np.array(casepermutations[protID])
                controlrecord=(record[0],record[1],controlname)
                overallresults=controlvalues[controlrecord]-self.clumps[record]
                diffdist=controls-cases
                pgreater=float(len(diffdist[diffdist<overallresults]))/float(len(diffdist))
                pless=float(len(diffdist[diffdist>overallresults]))/float(len(diffdist))
                print '\t'.join(record)+'\t'+str(overallresults)+'\t'+str(pgreater)+'\t'+str(pless)+'\t'+str(self.clumps[record])+'\t'+str(controlvalues[controlrecord])
            except KeyError:
                pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main():


    opts=CmdOpts()
    MyParser = Parser(opts.filename,opts.proteinlengths,float(opts.affreq),opts.nmutcutoff,opts.normalize)
    #MyParser()
    if(opts.controlperms and opts.permutations):
        if opts.titles:
            print 'GENE\tPROTEIN_ID\tSTUDY_NAME\tCLUMP_SCORE_DIFFERENCE(Controls-Cases)\tP-value(Probability Cases have a CLUMP score greater than the Controls)\tP-value(Probability Cases have a CLUMP score less than Controls)\tCASES_Raw_Clump_Score\tCONTROL_Raw_Clump_Score'
        MyParser.performpermutations(opts.filename,opts.controlperms,int(opts.permutations),float(opts.affreq))
    else:
        if opts.titles:
            print "GENE\tPROTEIN_ID\tSTUDY_NAME\tRaw_Clump_Score"
        MyParser()
if __name__ == '__main__':
    main()
