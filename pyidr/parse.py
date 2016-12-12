#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) <year> University of Dundee & Open Microscopy Environment.
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""
Module documentation
"""

#######################################################################
# create_bulk_annotations_file_using_study.pl
#
# Eleanor Williams 2015-06-08
#
# Script to read in study, library file and processed data
# files for a High Content Screen and output a comma separated file that can be
# used as a bulk annotation of the screen in omero.
#
#######################################################################

#######################################################################
# what it does
#
# A. reads in a study, library and processed data file (hit list)
# plus the number of the screen that you want to generate the bulk
# annotation file for e.g. 1, 2, 3.
# Usually screenA = 1, screenB = 2, screenC = 3 etc
#
# B. from the study file the screen get:
#   - the rows relating to the screen in question
#   - the column heading to be used to combine the library and
#     processed data files e.g. Gene Identifier
#   - the phenotypes (if any)
#   - the URIs that should be added to the bulk annotation column headings
#
# C. from the library file find out:
#   - which column contains the identifier that is going to be used to
#     combined the data with the the processed file
#
# D. from the library and processed data files:
#   - how many columns are in common between the library and processed file
#   - work out how many blank columns will need to be added if there is
#   processed data for a well
#
# E. from the processed file:
#    - find out which column contains the identifier which is going to
#    match up with the library file
#    - remove the columns from the processed file that are already in
#    the library file as don't need them twice in the final bulk annotation
#    file
#    - for each phenotype in the processed file find out what the
#    associated ontology mappings are (if any)
#
# F. print out data to an output file.
#    - goes through each line of the library file
#    - and adds in the processed data with the ontology mappings added
#      (if there is any)
#    - if there is no processed data then adds blank columns
#    - prints out each line to an output file
#
#######################################################################

#######################################################################
# Things to be careful of
#
# Line Endings
# Watch out for mac line endings in the files. Needs to be unix line
# endings.
#

#######################################################################

#######################################################################
# TODOs
#
# - if join on Plate_Well then remove this column from the bulk annoation
#   as don't really need it there.
# - put a lot more into subroutines esp the phenotype ontology mapping
#   part.
# - could do a lot more checking of values e.g. check URIs and
#   ontology accessions are correct format etc.
# - probably lots of things could be simplified and improved on
#
#######################################################################

#######################################################################
# REVISIONS
#
# 11-12-2016
#
# 1. Rewrite in Python
#
# 14-06-2016
#
# 1. removed adding URL into column headings


######################################################################
# A. get inputs from user and open the files
######################################################################

from re import match
from argparse import ArgumentParser
from logging import basicConfig
from logging import WARN
from logging import getLogger

LOG = getLogger("idr.parse")


def cli():
    parser = ArgumentParser(description=(
        "Creates a bulk annotation file for a HCS in Omero from a library "
        "file and processed data file.\n\n"
        "Example:  %s idr0000-study.txt idr0000-screenB-library.txt "
        "idr0000-screenB-processed.txt 2\n"
        "Output :  The output file name is taken from the library file with "
        "the extension -annotation.txt rather than -library.txt"
    ))
    parser.add_argument("-v", action="count", default=0,
                        help="increase verbosity")
    parser.add_argument("studyFile")
    parser.add_argument("libraryFile")
    parser.add_argument("processedDataFile",
                        help="hit list")
    parser.add_argument("screenNumber", type=int,
                        help="(1,2,3 etc)")
    # TODO: likely this should just be a idr directory and screen number
    ns = parser.parse_args()
    lvl = WARN - (ns.v * 10)
    basicConfig(level=lvl, format='%(message)s')
    sf = StudyFiles(ns.studyFile, ns.libraryFile, ns.processedDataFile)
    sf.read_files()
    col, phenos = sf.process_screen(ns.screenNumber)
    sf.process_library(col, phenos)
    return sf  # Returned for testing


class StudyFiles(object):

    def __init__(self, study, library, processed_data):
        self.study = study
        self.library = library
        self.processed_data = processed_data

    def read_files(self):
        """
        Read in the entirety of the configured files.
        """
        with open(self.study) as o:
            self.study_contents = o.read()
        with open(self.library) as o:
            self.library_lines = o.readlines()
        with open(self.processed_data) as o:
            self.processed_lines = o.readlines()
        # TODO: these should be generators


######################################################################
# B. process the study file
# 1. find the right screen
# 2. find which column to combine on
# 3. get the phenotypes and the ontology mappings
#####################################################################

    def process_screen(self, screenNumber):

        # 1. find the right screen
        self._parseScreenSections(screenNumber)
        section = self.screen_sections[screenNumber-1]
        for row in section[0:5]:
            LOG.debug(">>>: %s" % row)
        screenName = section[1].split("\t")[1]
        LOG.info("Making bulk annotation file for screen %s" % screenName)

        # 2. find which column to combine on
        columnTitleToCombineOn = self._getColumnToCombineOn(section)
        LOG.info("The column to combine on is %s" % columnTitleToCombineOn)

        # 3. get the phenotypes and the ontology mappings
        # Create a hash of the submitter phenotypes and their ontology sources,
        # terms and accessions
        # There may be more than one ontology term so structures will be like
        # submitted_phenotype1 => CMPO, CMPO_term, CMPO_acccession
        # submitted_phenotype2 => CMPO, CMPO_term, CMPO_acccession, CMPO,
        # CMPO_term, CMPO_acccession

        phenotype_ontologyArray = self._getPhenotypes(section)
        if not phenotype_ontologyArray:
            LOG.warn("WARNING: There are no phenotypes for this screen")

        # 4. now ommitted (was getting URLS)

        return columnTitleToCombineOn, phenotype_ontologyArray

######################################################################
# C. process the library file
# 5. find out which column contains the identifier which is going to
#    match up with the processed file
######################################################################

    def process_library(self, columnTitleToCombineOn, phenotype_ontologyArray):

        # 5. which column in library file has identifier to match column in
        # processed file

        indexOfLibraryFileColumnForMatching = None

        libraryHeaderRow = [x.strip()
                            for x in self.library_lines[0].split("\t")]

        for n, column in enumerate(libraryHeaderRow):
            # TODO: clarify regex
            if column == columnTitleToCombineOn:
                indexOfLibraryFileColumnForMatching = n
                break

        LOG.info("index of library file column for matching is %s" %
                 indexOfLibraryFileColumnForMatching)

        # remove any new line characters
        # TODO (perhaps already done)

####################################################################
# D. columns in common between library and processed data files
# (note no effort is made to check the content is the same, just
# the column titles)
#
# 6. find out which columns appear in both the library and processed
#    files? We don't want them repeated twice in the final output
#    file.
# 7. if a well in the library file has no processed data information
#    then we will need to add a number of empty columns to that row
#    in the output file so that table has same number of columns for
#    each row.  So work out how many blank columns to add.
######################################################################

        # 6. which columns appear in both the library and processed files?

        columnsToLooseFromProcessedFile = set()
        numberOfColumnsUniqueToProcessedFile = 0
        blankColumnsIfNoProcessedData = ""

        processedHeaderRow = self.processed_lines[0].split("\t")
        for index, col in enumerate(processedHeaderRow):
            col = col.strip()
            LOG.debug("Column is %s" % processedHeaderRow[index])
            # TODO: review regex
            #       have to do quotemeta to match if the string has square
            #       brackets  e.g. Experimental Condition [genotype]
            if col in libraryHeaderRow:
                columnsToLooseFromProcessedFile.add(index)
        LOG.debug("columnsToLooseFromProcessedFile: %s" %
                  columnsToLooseFromProcessedFile)

        numberOfColumnsUniqueToProcessedFile = len(processedHeaderRow) - \
            len(columnsToLooseFromProcessedFile)

        # 7. make a string of blank columns equal to the number of columns
        #    left in the processed data file after removing columns also in
        #    the library file

        for blanks in range(numberOfColumnsUniqueToProcessedFile-1):
            blankColumnsIfNoProcessedData += ","

######################################################################
# E. process the processed data file
# 8. find out which column contains the identifier which is going to
#    match up with the library file
# 9. remove the columns from the processed file that are already in
#    the library file as don't need them twice in the final file
# 10. For each phenotype in the processed file find out what the
#    associated ontology mappings are (if any)
######################################################################

        # 8. which column in processed file has identifier to match column in
        # library file

        indexOfProcessedFileColumnForMatching = None

        for p, column in enumerate(processedHeaderRow):
            # have to quotemeta this as column to match
            # on might have brackets in it
            # e.g. Experimental Condition [cell line]
            if column == columnTitleToCombineOn:
                LOG.debug("column to match on is %s (%s)" % (column, p))
                indexOfProcessedFileColumnForMatching = p
                break

        # 9. remove the columns from the processed file that are already in
        #    the library file as don't need them twice in the final file. To
        #    do this need to create a hash of the column numbers and the
        #    values otherwise as soon as one column is removed all the other
        #    column numbers in the array will change.
        #    Then put the remaining row into a hash with the common
        #    identifier as the key.

        Identifier_otherColumns = dict()

        # get each row of the processed file
        for row in self.processed_lines:
            row = row.strip()
            row = row.split("\t")

            # then create new array with just the column values we want to keep
            thisRowColumnValuesToKeep = []
            for key, val in enumerate(row):
                if key not in columnsToLooseFromProcessedFile:
                    thisRowColumnValuesToKeep.append(val)

            idx = row[indexOfProcessedFileColumnForMatching]
            Identifier_otherColumns[idx] = thisRowColumnValuesToKeep

        # 10. For each phenotype in the processed file find out what the
        #     associated ontology mappings are (if any)

        # first clone the hash of arrays with the identifier and the columns so
        # can move through all the columns in the original file without the
        # column numbers jumping due to column insertions

        Identifier_otherColumnsWithOntology = dict(Identifier_otherColumns)

        LOG.info("At start old header row is %s" %
                 Identifier_otherColumns[columnTitleToCombineOn])
        LOG.info("At start new header row is %s" %
                 Identifier_otherColumnsWithOntology[columnTitleToCombineOn])

        # if there are any phenotypes then add the ontology to the new header
        # row otherwise it will just stay the same as it is

        # if there are any phenotypes mentioned in the study file for the screen
        if phenotype_ontologyArray:

            # find which are the phenotype columns from the processed file
            # when get a phenotype column, find out what its in it, then get
            # mapping store which ontology with the phenotype to add link in
            # header then add the header rows then go through and add the
            # mappings

            # going through column headings of original array
            sz = len(Identifier_otherColumns[columnTitleToCombineOn])
            for a in range(sz):
                # Reminder: %Identifier_otherColumns has each identifier e.g.
                # Plate_Well or Gene Identfier that is used to combined the
                # library and processed data files as the key, and all the
                # processed data columns that go with that identifier as the
                # values So here we are going through the column headings row
                # because the key is what ever the column title to combine on is

                # when we get a phenotype column ...
                val = Identifier_otherColumns[columnTitleToCombineOn][a]
                if match("^Phenotype\s?\d*", val):

                    mapping = list()
                    ontologiesUsed = list()

                    # FIRST TIME ROUND - JUST FIND OUT FOR THIS PHENOTYPE IF
                    # THERE IS A MAPPING, AND IF SO IS IT ONE OR TWO TERMS
                    # AND WHAT ONTOLOGIES ARE THEY FROM

                    for identifier in Identifier_otherColumns.keys():
                        # start going through all the rows in that phenotype
                        # column to find one with a value

                        LOG.error("identifier: %s" % identifier)
                        if (Identifier_otherColumns[identifier][a]
                                and identifier != columnTitleToCombineOn):

                            # if the value matches a word character but is not
                            # the column heading; see if this phenotype has an
                            # ontology mapping

                            # check the phenotype exists in the study file.
                            # Have to use quotemeta in case there is a bracket
                            # or other character needing escaping in the
                            # phenotype value (TODO)
                            if Identifier_otherColumns[identifier][a] in \
                                    phenotype_ontologyArray.keys():

                                # this info comes from the study file
                                mapping = phenotype_ontologyArray[
                                    Identifier_otherColumns[identifier][a]]

                                if len(mapping) == 3:
                                    ontologiesUsed.append(mapping[0])
                                    break
                                elif len(mapping) == 6:
                                    ontologiesUsed.append(mapping[0])
                                    ontologiesUsed.append(mapping[3])
                                    break
                            else:
                                raise Exception((
                                    "ERROR: Phenotype '%s' does not exist in "
                                    "the study file") %
                                    Identifier_otherColumns[identifier][a]
                                )

        # get the URIs associated with the ontologies used (REMOVED)

        """
         # SECOND TIME ROUND - GO THROUGH EACH IDENTIFIER, IF THERE IS A PHENOTYPE ADD THE MAPPINGS, IF NOT THEN JUST ADD THE SAME NUMBER OF TABS
     # IF THERE IS A PHENOTYPE BUT THIS PHENOTYPE HAS NO ONTOLOGY MAPPING THEN SKIP TO THE NEXT COLUMN

         if ($numberOfMappings == 1){

              # for each identifier, if there is a value - put in the one mapping, if heading, put in column titles, if no value put in 2 tabs
       foreach my $identifier (keys %Identifier_otherColumns){

                   if (($Identifier_otherColumns{$identifier}[$a] =~ m/\w+/) && ($identifier ne $columnTitleToCombineOn))  { # if the value matches a word character but is not the column heading

                 @mapping = @{$phenotype_ontologyArray{$Identifier_otherColumns{$identifier}[$a]}};
                     splice @{$Identifier_otherColumnsWithOntology{$identifier}}, $b+1, 0, $mapping[1], $mapping[2];

               }elsif($identifier eq $columnTitleToCombineOn){ # column title
                           if(${$Identifier_otherColumns{$columnTitleToCombineOn}}[$a] =~ m/^Phenotype\s?(\d+)$/){# if column name has a number e.g. Phenotype 1
                 my $number = $1;
                 my $termName = "Phenotype ".$number." Term Name";
                 my $termAcc = "Phenotype ".$number." Term Accession";
                # my $termAcc = "Phenotype ".$number." Term Accession %% url=".$ontology_URI{$ontologiesUsed[0]}."%s";

                 splice @{$Identifier_otherColumnsWithOntology{$columnTitleToCombineOn}}, $b+1, 0, $termName, $termAcc;
                 $numberOntologyColumnsAdded = $numberOntologyColumnsAdded + 2;
                         }else{
                 my $termAcc = "Phenotype Term Accession";
                 #my $termAcc =  "Phenotype Term Accession %% url=".$ontology_URI{$ontologiesUsed[0]}."%s";
                           splice @{$Identifier_otherColumnsWithOntology{$columnTitleToCombineOn}}, $b+1, 0, 'Phenotype Term Name', $termAcc;
                       }

             }else{ # no value
             # insert empty column
                     splice @{$Identifier_otherColumnsWithOntology{$identifier}}, $b+1, 0,"", "";

                    }
               } # foreach identifier
         $b=$b+3;
     }elsif($numberOfMappings == 2){
                # for each identifier, if there is a value - put in the two mappings, if heading, put in column titles, if no value put in 4 tabs

          foreach my $identifier (keys %Identifier_otherColumns){
        if (($Identifier_otherColumns{$identifier}[$a] =~ m/\w+/) && ($identifier ne $columnTitleToCombineOn))  { # if the value matches a word character but is not the column heading
                @mapping = @{$phenotype_ontologyArray{$Identifier_otherColumns{$identifier}[$a]}};
                    splice @{$Identifier_otherColumnsWithOntology{$identifier}}, $b+1, 0, $mapping[1], $mapping[2], $mapping[4], $mapping[5];

         }elsif($identifier eq $columnTitleToCombineOn){ # column title
                          if(${$Identifier_otherColumns{$columnTitleToCombineOn}}[$a] =~ m/^Phenotype\s?(\d+)$/){# if column name has a number e.g. Phenotype 1
                     my $number = $1;
                 my $termNameA = "Phenotype ".$number." Term Name a";
                             my $termAccA = "Phenotype ".$number." Term Accession a";
                 #my $termAccA = "Phenotype ".$number." Term Accession a %% url=".$ontology_URI{$ontologiesUsed[0]}."%s";
                          my $termNameB = "Phenotype ".$number." Term Name b";
                 my $termAccB = "Phenotype ".$number." Term Accession b";
                 #my $termAccB = "Phenotype ".$number." Term Accession b %% url=".$ontology_URI{$ontologiesUsed[1]}."%s";
                             splice @{$Identifier_otherColumnsWithOntology{$columnTitleToCombineOn}}, $b+1, 0, $termNameA, $termAccA, $termNameB, $termAccB;
               }else{
                 my $termAccA = "Phenotype Term Accession a";
                             my $termAccB = "Phenotype Term Accession b";
                             #my $termAccA = "Phenotype Term Accession a %% url=".$ontology_URI{$ontologiesUsed[0]}."%s";
                             #my $termAccB = "Phenotype Term Accession b %% url=".$ontology_URI{$ontologiesUsed[1]}."%s";
                             splice @{$Identifier_otherColumnsWithOntology{$columnTitleToCombineOn}}, $b+1, 0, 'Phenotype Term Name a', $termAccA, 'Phenotype Term Name b', $termAccB;
                      }
                 $numberOntologyColumnsAdded = $numberOntologyColumnsAdded + 4;
         }else{ # no value
                         # inserting four spaces
                     splice @{$Identifier_otherColumnsWithOntology{$identifier}}, $b+1, 0,"", "", "", "";
                    }

               } # foreach identifier
         $b=$b+5;
     }else{
       # must be no mappings for the phenotype so do nothing
     $b++;
     }


       }else{ # not a phenotype column so just move to the next column
     $b++;
       }

} # for each column heading




#print "At end old header row is @{$Identifier_otherColumns{$columnTitleToCombineOn}}\n";
#print "At end new header row is @{$Identifier_otherColumnsWithOntology{$columnTitleToCombineOn}}\n";
#print "Total number of columns added due to ontologies is $numberOntologyColumnsAdded\n";

# add on the number of columns added for ontologies to the list of columns to be added if there is no processed data

for (my $x=0; $x<$numberOntologyColumnsAdded; $x++){
$blankColumnsIfNoProcessedData = $blankColumnsIfNoProcessedData.",";
}

} # if there are any phenotypes listed in the study file

######################################################################
# F. create the output file                                             #
# 11. create the out file name and open it                            #
# 12. the library file contains information for every well, so go     #
#    through each line in it, print out each line, adding            #
#    information from the processed file if there is any.            #
######################################################################

######################################################################
# 11. open the output file

my $outfile = $libraryFile;
$outfile =~ s/-library\.txt/-annotation\.txt/g; # its comma delimited but
#                                                 make it end in .txt so can
#                                                 open in Excel corrently
open (OUT, ">$outfile");

######################################################################
# 12. go through each line in the library file and print out          #
#    adding processed file info if needed                            #

my $v=0;
foreach my $libRow (@libraryFile){
  chomp($libRow);

  $libRow =~ s/\t/\,/g; # change tabs to commas as want output file to be comma separated
  print OUT "$libRow\,";

  my @libraryRow = split ("\,", $libRow, -1);



 # add processed data if there is any or blank columns if not        #

    if (exists ($Identifier_otherColumnsWithOntology{$libraryRow[$indexOfLibraryFileColumnForMatching]})) {
      my $processedRow = join("\,", @{$Identifier_otherColumnsWithOntology{$libraryRow[$indexOfLibraryFileColumnForMatching]}});
      print OUT "$processedRow";
    }else{
      print OUT "$blankColumnsIfNoProcessedData";
    }


   # put a line ending at the end of every row
   print OUT "\n";

  $v++;
}

close (OUT);
"""

######################################################################
# Helpers
######################################################################

    def _parseScreenSections(self, screenNumber):
        if "Screen Number" not in self.study_contents:
            raise Exception(("Phrase 'Screen Number' does not exist in the "
                             "study file so can't get the screen information: "
                             "%s") % self.study)

        # split the studyFile on "Screen Number"
        sections = self.study_contents.split("Screen Number")
        self.study_section = sections[0].split("\n")
        self.screen_sections = [s.split("\n") for s in sections[1:]]

        # first section [0] is the study top level info then each following
        # section is a screen. get the rows for the screen we want as long as
        # the screen number is valid
        found = len(self.screen_sections)
        if found < screenNumber:
            raise Exception((
                "Screen Number must be the same as one of those specified in "
                "the study file: %s > %s") % (screenNumber, found))

    def _getColumnToCombineOn(self, section):
        # Find out which column the library and processed files should be
        # joined on from the study file, screen section

        columnToCombine = None
        for row in section:
            if "Processed Data Column Link To Library File" in row:
                cells = row.split("\t")
                # the column to combine on is always listed in the first column
                # after the 'Processed Data Column Link To Library File' tag
                columnToCombine = cells[1]
                break

        # if there is no 'Processed Data Column Link to Library File' row
        # or the value is empty then stop here
        if not columnToCombine:
            raise Exception((
                "No column to combine on information for the screen "
                "in the study file"))

        return columnToCombine

    def _getPhenotypes(self, section):

        screenPhenotype_ontologyArray = {}
        # keep it easy to understand using variables
        phenotypes = []
        termSource1 = []
        termNames1 = []
        termAccs1 = []
        termSource2 = []
        termNames2 = []
        termAccs2 = []

        # find where the phenotype row is
        # then get the first mappings if there are any
        # then get the second mappings if there are any
        # assuming there are never more than 2 mappings

        for n, row in enumerate(section):

            if "Phenotype Name" in row:
                # split the row and get all the phenotypes
                # + first lot of mappings (always have these rows)
                phenotypes = row.split("\t")
                termSource1 = section[n+3].split("\t")
                # TODO: check the -1 keeps trailing columns
                termNames1 = section[n+4].split("\t")
                termAccs1 = section[n+5].split("\t")

                # then see if we have another lot of ontology terms
                if "Phenotype Term Source REF" in section[n+6]:
                    LOG.debug("we have some second ontology terms")
                    termSource2 = section[n+6].split("\t")
                    termNames2 = section[n+7].split("\t")
                    termAccs2 = section[n+8].split("\t")

        LOG.debug("Getting the phenotypes from the library file")
        # put it all into the hash of arrays
        # maxNumberOfOntTerms = 1 (unused)
        # use this to see if we have any phenotypes
        # that map to 2 ontology terms
        # like 'abnormal' + 'microtubule structure'
        # go through each column in each row
        # first column is 'Phenotype Name' so don't need that
        for t in range(1, len(phenotypes)):
            # empty colection of items in one column
            mappingsArray = []
            LOG.debug("phenotype found is %s" % phenotypes[t])
            # only add values if they have a letter character i.e. not blank
            for arr in (termSource1, termNames1, termAccs1):
                if arr[t].strip():
                    mappingsArray.append(arr[t])

            if termSource2:
                LOG.debug("Looking at second mappings")
                for arr in (termSource2, termNames2, termAccs2):
                    if arr[t].strip():
                        mappingsArray.append(arr[t])

            LOG.debug("array is going to be %s" % mappingsArray)
            screenPhenotype_ontologyArray[phenotypes[t]] = mappingsArray

        return screenPhenotype_ontologyArray


if __name__ == "__main__":
    sf = cli()
