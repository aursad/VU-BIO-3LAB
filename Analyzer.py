__author__ = 'Aurimas Sadauskas'
_version__ = '1.0'
_lastUpdate = '2014-12-15'

import os

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Timer import Timer
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class Analyzer:
    # Constructor
    def __init__(self, program, database):
        self.program = program
        self.record = ""
        self.database = database
        self.cdHittedFilesArray = []

    def startBlast(self, entrezQuery, input="search_output.xml"):
        blastRecords = self.blast(entrezQuery)
        fileToWrite = open(input, "w")
        arrayOfSequences = []
        for sequence in blastRecords.alignments:
            #changing seq length be more than 60
            if len(sequence.hsps[0].sbjct) > 60:
                sequenceToAppend = SeqRecord(Seq(sequence.hsps[0].sbjct), id=sequence.title)
                arrayOfSequences.append(sequenceToAppend)
        SeqIO.write(arrayOfSequences, fileToWrite, "fasta")
        fileToWrite.close()

        return blastRecords;

    def blast(self, entrezQuery):
        blastOutput = NCBIWWW.qblast(self.program, self.database, self.record.seq,
                                     entrez_query=entrezQuery, hitlist_size=1000, expect=100.0)
        return NCBIXML.read(blastOutput);

    def loadMainSeq(self, path, format="fasta"):
        #loading main genom fragment sequence
        self.record = SeqIO.read(open(path), format=format)

    def startMafft(self, input, output="mafft_output.fasta"):
        os.system("mafft --quiet {0} > {1}".format(input, output))
        return;

    def startCdHit(self, inputBlast, inputCdHit):
        self.cdHittedFilesArray += [inputCdHit]
        #running cdhit from specific file destination
        os.system('cd-hit-windows\cd_hit -i {0} -o {1}'.format(inputBlast, inputCdHit))

    def joinFiles(self, tempFile):
        fileForMerge = open(tempFile, "w")
        #Joining all files to one, from property for cd hitted files array list
        for fileName in self.cdHittedFilesArray:
            fileForMerge.write(open(fileName).read())
        fileForMerge.close()

def main():
    #List of dangerious types
    dangerious = [16, 18, 31, 33, 35, 51, 52]
    #List of not dangerious types
    notDangerious = [6, 11, 40, 42, 43, 44, 57, 81]
    #Final file after mafft
    mafftOutputFile = "mafft_output.fasta"
    #File output for temp joined cd hitted files
    mafftTempFile = "mafft_temp.fasta"
    #Directory where is saving all blast search files
    blastFolder = "blast"
    #Directory where is saving all file after cdhit
    cdHitFolder = "cdhit"
    #Prefix of blast files
    fileNameForSequences = "blast_type_"
    #Prefix of cdHit files
    fileNameForSequencesAfterCdHit = 'cdhit_type_'

    try:
        with Timer() as timer:
            #Check if files dir's exist if not to create
            if not os.path.exists(blastFolder):
                os.makedirs(blastFolder)
            if not os.path.exists(cdHitFolder):
                os.makedirs(cdHitFolder)

            #Start Seq analyzer with constructor
            seqAn = Analyzer("blastn", "nr")
            #Loading main genom fragment
            seqAn.loadMainSeq("human_papillomavirus_HPV16.fasta")

            #Generating list of types, first dangerious after then not dangerious getting
            for virus in dangerious + notDangerious:
                #Generating blast search query
                enterezQuery = '"papillomavirus"[Organism] AND ( *"type %s"[title] AND *human*[title])' % virus

                #Creating full path to save blast result
                completeNameBlast = os.path.abspath("{0}/{1}{2}.fa".format(blastFolder, fileNameForSequences, virus.__str__()))
                print('Processing file: {0}'.format(completeNameBlast.__str__()))

                #Doing blast for specific query and then file is saved
                if not os.path.isfile(completeNameBlast):
                    seqAn.startBlast(enterezQuery, completeNameBlast)

                #Creating full path to save cdHit file of specific virus
                completeNameCdHit = os.path.abspath("{0}/{1}{2}.fa".format(cdHitFolder, fileNameForSequencesAfterCdHit, virus.__str__()))

                #if cdhit file for type not exist doing cdhit else file name saved to list of files
                if not os.path.isfile(completeNameCdHit):
                    seqAn.startCdHit(completeNameBlast, completeNameCdHit)
                else:
                    seqAn.cdHittedFilesArray += [completeNameCdHit]

                print('Processing file {0} end.'.format(completeNameBlast.__str__()))

            # Remove joined file if exist
            if os.path.isfile(mafftTempFile):
                os.remove(mafftTempFile)
            # Remove mafft output file if exist
            if os.path.isfile(mafftOutputFile):
                os.remove(mafftOutputFile)

            # Joining all cdhit files and then starting mafft
            seqAn.joinFiles(mafftTempFile)
            #start mafft for joined file, saved to mafftOutputFile
            seqAn.startMafft(mafftTempFile, mafftOutputFile)
    finally:
        #Showing time for how long program runned
        print('Request took %.03f sec.' % timer.interval)

# Run main
main()