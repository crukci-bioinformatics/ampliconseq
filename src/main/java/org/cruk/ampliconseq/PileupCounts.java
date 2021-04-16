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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;
import org.cruk.htsjdk.intervals.IntervalUtils;
import org.cruk.htsjdk.pileup.PileupUtils;

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

/**
 * Command line tool for generating a pileup summary for all loci within a
 * specified set of intervals.
 *
 * Each interval is processed separately so that the interval name can be
 * included in the output.
 *
 * @author eldrid01
 */
public class PileupCounts extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    private String id;
    private File bamFile;
    private File intervalsFile;
    private File referenceSequenceFile;
    private File pileupCountsFile;
    private int minimumBaseQuality = 10;
    private int minimumMappingQuality = 1;

    public static void main(String[] args) {
        PileupCounts pileupCounts = new PileupCounts();
        pileupCounts.parseCommandLineArgs(args);
        pileupCounts.run();
    }

    @Override
    protected String getHelpDescription() {
        return "Generates a pileup summary with read counts for each position and allele";
    }

    @Override
    protected Options createOptions() {
        Options options = super.createOptions();

        Option option = new Option(null, "id", true,
                "Identifier for this dataset; if included the pileup counts table will have an additional ID column (required)");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("i", "input", true,
                "Input BAM file which must be in coordinate sort order and indexed (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("l", "intervals", true,
                "Intervals over which to generate pileup counts; can be in BED or Picard-style interval format (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("r", "reference-sequence", true, "Reference sequence FASTA file (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("o", "output", true,
                "The output pileup summary table with read counts for each position and allele (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("q", "minimum-base-quality", true,
                "Minimum base quality for the base call at a locus for reads to be included (default: "
                        + minimumBaseQuality + ")");
        option.setType(Number.class);
        options.addOption(option);

        option = new Option("m", "minimum-mapping-quality", true,
                "Minimum mapping quality for reads to be included (default: " + minimumMappingQuality + ")");
        option.setType(Number.class);
        options.addOption(option);

        return options;
    }

    @Override
    protected void extractOptionValues(CommandLine commandLine) throws ParseException {
        id = commandLine.getOptionValue("id");
        bamFile = (File) commandLine.getParsedOptionValue("input");
        intervalsFile = (File) commandLine.getParsedOptionValue("intervals");
        referenceSequenceFile = (File) commandLine.getParsedOptionValue("reference-sequence");
        pileupCountsFile = (File) commandLine.getParsedOptionValue("output");
        if (commandLine.hasOption("minimum-base-quality")) {
            minimumBaseQuality = ((Number) commandLine.getParsedOptionValue("minimum-base-quality")).intValue();
        }
        if (commandLine.hasOption("minimum-mapping-quality")) {
            minimumMappingQuality = ((Number) commandLine.getParsedOptionValue("minimum-mapping-quality")).intValue();
        }
    }

    /**
     * Main run method for iterating over loci within the given intervals and
     * tabulating read counts for each base.
     */
    private void run() {

        ProgressLogger progress = new ProgressLogger(logger, 100, "loci");

        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(intervalsFile);
        IOUtil.assertFileIsReadable(referenceSequenceFile);
        IOUtil.assertFileIsWritable(pileupCountsFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(bamFile);

        ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(
                referenceSequenceFile);

        List<Interval> intervals = IntervalUtils.readIntervalFile(intervalsFile);

        BufferedWriter writer = IOUtil.openFileForBufferedWriting(pileupCountsFile);

        try {
            writeHeader(writer);

            for (Interval interval : intervals) {

                logger.info("Interval: " + interval.getName());

                IntervalList intervalList = new IntervalList(referenceSequenceFileWalker.getSequenceDictionary());
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

            CloserUtil.close(reader);

            writer.close();

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
