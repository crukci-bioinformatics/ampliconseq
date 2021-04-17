/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.ampliconseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;
import org.cruk.htsjdk.intervals.IntervalUtils;
import org.cruk.htsjdk.pileup.PileupUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Command line tool for generating a pileup summary for all loci within a
 * specified set of intervals.
 *
 * Each interval is processed separately so that the interval name can be
 * included in the output.
 *
 * @author eldrid01
 */
@Command(name = "pileup-counts", versionProvider = PileupCounts.class, description = "\nGenerates a pileup summary with read counts for each position and allele.\n", mixinStandardHelpOptions = true)
public class PileupCounts extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    @Option(names = "--id", required = true, description = "Identifier for this dataset included in the output table in the ID column (required).")
    private String id;

    @Option(names = { "-i",
            "--input" }, required = true, description = "Input BAM file which must be in coordinate sort order and indexed (required).")
    private File bamFile;

    @Option(names = { "-l",
            "--intervals" }, required = true, description = "Intervals over which to generate pileup counts; can be in BED or Picard-style interval format (required).")
    private File intervalsFile;

    @Option(names = { "-r",
            "--reference-sequence" }, required = true, description = "Reference sequence FASTA file which must be indexed and have an accompanying dictionary (required).")
    private File referenceSequenceFile;

    @Option(names = { "-o",
            "--output" }, required = true, description = "The output pileup summary table with read counts for each position and allele (required).")
    private File pileupCountsFile;

    @Option(names = { "-q",
            "--minimum-base-quality" }, description = "Minimum base quality for the base call at a locus for reads to be included (default: ${DEFAULT-VALUE}).")
    private int minimumBaseQuality = 10;

    @Option(names = { "-m",
            "--minimum-mapping-quality" }, description = "Minimum mapping quality for reads to be included (default: ${DEFAULT-VALUE}).")
    private int minimumMappingQuality = 1;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new PileupCounts()).execute(args);
        System.exit(exitCode);
    }

    /**
     * Main run method for iterating over loci within the given intervals and
     * tabulating read counts for each base.
     */
    @Override
    public void run() {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        ProgressLogger progress = new ProgressLogger(logger, 100, "loci");

        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(intervalsFile);
        IOUtil.assertFileIsReadable(referenceSequenceFile);
        IOUtil.assertFileIsWritable(pileupCountsFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(bamFile);

        ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(
                referenceSequenceFile);

        SAMSequenceDictionary sequenceDictionary = referenceSequenceFileWalker.getSequenceDictionary();
        if (sequenceDictionary == null) {
            logger.error("Sequence dictionary for reference sequence is missing");
            System.exit(1);
        }

        List<Interval> intervals = IntervalUtils.readIntervalFile(intervalsFile);

        BufferedWriter writer = IOUtil.openFileForBufferedWriting(pileupCountsFile);

        try {
            writeHeader(writer);

            for (Interval interval : intervals) {

                logger.info("Interval: " + interval.getName());

                IntervalList intervalList = new IntervalList(sequenceDictionary);
                intervalList.add(interval);

                SamLocusIterator locusIterator = new SamLocusIterator(reader, intervalList);

                // exclude reads that are marked as failing platform/vendor quality checks
                locusIterator.setIncludeNonPfReads(false);

                // exclude secondary alignments and reads marked as duplicates
                locusIterator.setSamFilters(Arrays.asList(new SecondaryAlignmentFilter(), new DuplicateReadFilter()));

                SamLocusAndReferenceIterator locusAndReferenceIterator = new SamLocusAndReferenceIterator(
                        referenceSequenceFileWalker, locusIterator);

                for (SamLocusAndReferenceIterator.SAMLocusAndReference locusAndReference : locusAndReferenceIterator) {

                    List<RecordAndOffset> pileup = locusAndReference.getRecordAndOffsets();

                    List<RecordAndOffset> filteredPileup = PileupUtils.filterLowQualityScores(pileup,
                            minimumBaseQuality, minimumMappingQuality);

                    filteredPileup = PileupUtils.filterOverlappingFragments(filteredPileup);

                    writePileupCounts(writer, interval, locusAndReference, filteredPileup);

                    LocusInfo locusInfo = locusAndReference.getLocus();
                    progress.record(locusInfo.getContig(), locusInfo.getPosition());
                }

                locusAndReferenceIterator.close();
                locusIterator.close();
            }

            writer.close();

            referenceSequenceFileWalker.close();
            CloserUtil.close(reader);

        } catch (IOException e) {
            logger.error(e);
            System.exit(1);
        }

        logger.info("Finished");
    }

    /**
     * Write the header for the coverage/pileup output table.
     *
     * @param writer
     * @throws IOException
     */
    private void writeHeader(BufferedWriter writer) throws IOException {
        writer.write("ID");
        writer.write("\t");
        writer.write("Interval");
        writer.write("\t");
        writer.write("Chromosome");
        writer.write("\t");
        writer.write("Position");
        writer.write("\t");
        writer.write("Reference base");
        writer.write("\t");
        writer.write("Depth unfiltered");
        writer.write("\t");
        writer.write("Depth");
        writer.write("\t");
        writer.write("A count");
        writer.write("\t");
        writer.write("C count");
        writer.write("\t");
        writer.write("G count");
        writer.write("\t");
        writer.write("T count");
        writer.write("\t");
        writer.write("N count");
        writer.write("\n");
    }

    private void writePileupCounts(BufferedWriter writer, Interval interval,
            SamLocusAndReferenceIterator.SAMLocusAndReference locusAndReference, List<RecordAndOffset> filteredPileup)
            throws IOException {

        writer.write(id);

        writer.write("\t");
        writer.write(interval.getName());

        LocusInfo locusInfo = locusAndReference.getLocus();

        writer.write("\t");
        writer.write(locusInfo.getContig());
        writer.write("\t");
        writer.write(Integer.toString(locusInfo.getPosition()));
        writer.write("\t");
        writer.write(locusAndReference.getReferenceBase());

        writer.write("\t");
        writer.write(Integer.toString(locusAndReference.getRecordAndOffsets().size()));

        writer.write("\t");
        writer.write(Integer.toString(filteredPileup.size()));

        int[] baseCounts = PileupUtils.getBaseCounts(filteredPileup);
        for (int i = 0; i < baseCounts.length; i++) {
            writer.write("\t");
            writer.write(Integer.toString(baseCounts[i]));
        }

        writer.write("\n");
    }
}
