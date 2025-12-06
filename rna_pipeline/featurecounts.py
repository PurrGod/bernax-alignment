# featurecounts.py
# Group Members: Akhilesh Nidamanuri, David Jiricek, Evan Fitzhugh

"""
Function: Run featureCounts on all BAM files and parse the per-read
          assignment BAMs (-R BAM) into per-sample assignment summaries.
"""

#Imports important files including paths, dataclasses, tools relevant to featureCounts
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Dict, Any, Set
from .star_runner import StarBatchOutputs
from . import utils
import pysam

#This class is the named categories that you need for featureCounts
class AssignmentCategory (Enum):
    ASSIGNED = auto() #Read got assigned to a feature (gene).
    UNASSIGNED_UNMAPPED = auto() #Read wasn’t mapped at all.
    UNASSIGNED_NO_FEATURES = auto() #Read mapped, but not overlapping any annotated feature
    UNASSIGNED_MAPPING_QUALITY = auto() #The mapping quality too low.
    UNASSIGNED_AMBIGUITY = auto() #Read overlapped multiple genes/features

#This class corresponds to a collection of outputs from featureCounts
#Creates an object that contains all the outputs
@dataclass
class FeatureCountsResult :
    countsFile: Path
    summaryFile: Path
    perSampleAssignmentFiles: Dict[str, Path]

#This class represents what happens to a single sample
#The samples have assigned IDs, or unassigned IDs that are put into a dictionary
@dataclass
class SampleAssignments :
    assignedIds: set
    unassignedIds: set
    categoryCounts: Dict[AssignmentCategory, int]


def runFeatureCounts (starOutputs, args, refCfg, outDir):
    '''
    Function: runFeatureCounts
    Purpose: Run featureCounts on all BAM files.
    Outputs: FeatureCountsResult object
    '''
   #Ensures the directory exists
    utils.ensureDir(outDir)
    
    #Defines the files that it will create
    countsFile = outDir / "counts.txt"
    summaryFile = outDir / "counts.txt.summary"

    # Determine GTF from refCfg or args
    #If there are no GTF file it will raise an error
    gtf = refCfg.get("gtf") or getattr(args, "gtf", None)
    if gtf is None:
        raise ValueError(
            "No GTF file provided. Use --gtf or include 'gtf' in reference-config."
        )

    #Collect BAM files from STAR, raises an error if no BAM files
    bamFiles = [o.bam for o in starOutputs.perSample.values()]
    if not bamFiles:
        raise ValueError("No BAM files found in StarBatchOutputs.")

    #This just builds the commandline for featurecounts and runs it one time across the BAM files
    cmd = [
        "featureCounts",
        "-T", str(getattr(args, "threads", 4)),
        "-a", str(gtf),
        "-o", str(countsFile),
        "-R", "BAM",
        *(str(b) for b in bamFiles),
    ]
    #Prints the command, from utils
    utils.runCmd(cmd)

    # featureCounts with -R BAM creates <input>.featureCounts.bam
    # It saves them in the output directory (outDir), NOT the input directory
    perSampleAssignmentFiles = {}
    for sampleId, starOut in starOutputs.perSample.items():
        # The tool names the file: OriginalName.bam.featureCounts.bam
        outputBamName = starOut.bam.name + ".featureCounts.bam"
        perSampleAssignmentFiles[sampleId] = outDir / outputBamName

    #The bundle of results that it returns, in a single object
    return FeatureCountsResult(
        countsFile=countsFile,
        summaryFile=summaryFile,
        perSampleAssignmentFiles=perSampleAssignmentFiles,
    )

def parseAssignments (fcResult):
    '''
    Function: parseAssignments
    Purpose: Parses read assignment files into SampleAssignments.
    Outputs: Dictionary mapping sample ID to SampleAssignments
    '''
    
    #Where assignments will be stored
    assignments = {}

    #Basically this is the master organizing loop that opens each BAM files and checks it
    for sampleId, bamPath in fcResult.perSampleAssignmentFiles.items():
        
        #Creates the sets for categorizing the reads
        assignedIds = set()
        unassignedIds = set()
        categoryCounts = {
            cat: 0 for cat in AssignmentCategory
        }

        #Makes sure that the BAM file actually exists
        if not bamPath.exists():
            raise FileNotFoundError(f"Expected featureCounts output BAM not found at: {bamPath}")

        #For every BAM File reads it and checks the tag on it to categorize it
        with pysam.AlignmentFile(str(bamPath), "rb") as bam:
            for read in bam:
                readId = read.query_name
                tag = read.get_tag("XS") if read.has_tag("XS") else None
                
                #Categorizes the tags based on their labels and adds them to correct lists
                #Case 1 - assigned
                if tag == "Assigned":
                    assignedIds.add(readId)
                    categoryCounts[AssignmentCategory.ASSIGNED] += 1
                
                #Case 2 - unassigned
                else:
                    unassignedIds.add(readId)
                    if tag == "Unassigned_Unmapped":
                        categoryCounts[AssignmentCategory.UNASSIGNED_UNMAPPED] += 1
                    elif tag == "Unassigned_NoFeatures":
                        categoryCounts[AssignmentCategory.UNASSIGNED_NO_FEATURES] += 1
                    elif tag == "Unassigned_MappingQuality":
                        categoryCounts[AssignmentCategory.UNASSIGNED_MAPPING_QUALITY] += 1
                    elif tag == "Unassigned_Ambiguous":
                        categoryCounts[AssignmentCategory.UNASSIGNED_AMBIGUITY] += 1
        
        #Stores the data/results for each sample and moves on to the next
        assignments[sampleId] = SampleAssignments(
            assignedIds=assignedIds,
            unassignedIds=unassignedIds,
            categoryCounts=categoryCounts,
        )
    #Then returns the results
    return assignments