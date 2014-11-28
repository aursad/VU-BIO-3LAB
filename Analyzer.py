__author__ = 'Aurimas Sadauskas'
_version__ = '1.0'
_lastUpdate = '2014-11-27'

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
        self.entrezQuery = ""
        self.cdHittedFilesArray = []

    def startBlast(self, entrezQuery, input="search_output.xml"):
        entrezQuery = entrezQuery

        blastRecords = self.blast()
        fileToWrite = open(input, "w")
        arrayOfSequences = []
        for sequence in blastRecords.alignments:
            # padarom, kad pacios sekos ilgis butu didesnis uz 60
            if len(sequence.hsps[0].sbjct) > 60:
                sequenceToAppend = SeqRecord(Seq(sequence.hsps[0].sbjct), id=sequence.title)
                arrayOfSequences.append(sequenceToAppend)
        SeqIO.write(arrayOfSequences, fileToWrite, "fasta")
        fileToWrite.close()

        return blastRecords;

    def blast(self):
        blastOutput = NCBIWWW.qblast(self.program, self.database, self.record.seq,
                                     entrez_query=self.entrezQuery, hitlist_size=1000, expect=100.0)
        return NCBIXML.read(blastOutput);

    def loadMainSeq(self, path, format="fasta"):
        self.record = SeqIO.read(open(path), format=format)

    def startMafft(self, input, output="mafft_output.fasta"):
        os.system("mafft --quiet {0} > {1}".format(input, output))
        return;

    def startCdHit(self, inputBlast, inputCdHit):
        self.cdHittedFilesArray += [inputCdHit]
        os.system('cd-hit-windows\cd_hit -i {0} -o {1}'.format(inputBlast, inputCdHit))

    def joinFiles(self, tempFile):
        fileForMerge = open(tempFile, "w")
        for fileName in self.cdHittedFilesArray:
            fileForMerge.write(open(fileName).read())
        fileForMerge.close()

def main():
    dangerious = [16, 18, 31, 33, 35, 51, 52]
    notDangerious =  [6, 11, 40, 42, 43, 44, 57, 81]
    mafftOutputFile = "mafft_output.fasta"
    mafftTempFile = "mafft_temp.fasta"
    blastFolder = "blast"
    cdHitFolder = "cdhit"
    fileNameForSequences = "blast_type_"
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
            seqAn.loadMainSeq("human_papillomavirus_HPV16.fasta")

            for virus in dangerious + notDangerious:
                #Generating blast search query
                enterezQuery = '"papillomavirus"[Organism] AND ( *"type %s"[title] AND *human*[title])' % virus

                #Creating full path to save blast result
                completeNameBlast = os.path.abspath("{0}/{1}{2}.fa".format(blastFolder, fileNameForSequences, virus.__str__()))
                print('Processing file: {0}'.format(completeNameBlast.__str__()))

                if not os.path.isfile(completeNameBlast):
                    seqAn.startBlast(enterezQuery, completeNameBlast)

                #Creating full path to save cdHit file of specific virus
                completeNameCdHit = os.path.abspath("{0}/{1}{2}.fa".format(cdHitFolder, fileNameForSequencesAfterCdHit, virus.__str__()))
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
            seqAn.startMafft(mafftTempFile, mafftOutputFile)
    finally:
        print('Request took %.03f sec.' % timer.interval)

# Run main
main()