/**
 * MIT License
 *
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
 * Utility for generating a pileup summary for all loci within a specified set
 * of intervals.
 *
 * @author eldrid01
 */
public class PileupCounts extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    private String id;
    private File bamFile;
    private File ampliconsFile;
    private File referenceSequenceFile;
    private File pileupCountsFile;
    private int minimumBaseQuality = 0;
    private int minimumMappingQuality = 0;

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
                "Identifier for this dataset; if included the pileup counts table will have an additional ID column (optional");
        options.addOption(option);

        option = new Option("i", "bam", true,
                "BAM input file which must be in coordinate sort order and indexed (required)");
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

        option = new Option("o", "pileup-counts", true,
                "The output BAM file containing a pileup summary table with read counts for each position and allele (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("q", "minimum-base-quality", true,
                "Minimum base quality for the base call at a locus for reads to be included (default: "
                        + minimumBaseQuality + ")");
        option.setType(Number.class);
        options.addOption(option);

        option = new Option("m", "minimum-base-quality", true,
                "Minimum mapping quality for reads to be included (default: " + minimumMappingQuality + ")");
        option.setType(Number.class);
        options.addOption(option);

        return options;
    }

    @Override
    protected void extractOptionValues(CommandLine commandLine) throws ParseException {
        id = commandLine.getOptionValue("id");
        bamFile = (File) commandLine.getParsedOptionValue("bam");
        ampliconsFile = (File) commandLine.getParsedOptionValue("intervals");
        referenceSequenceFile = (File) commandLine.getParsedOptionValue("reference-sequence");
        pileupCountsFile = (File) commandLine.getParsedOptionValue("pileup-counts");
    }

    /**
     * Main run method for iterating over loci within the given intervals and
     * tabulating read counts for each base.
     */
    private void run() {

        ProgressLogger progress = new ProgressLogger(logger, 100, "loci");

        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(ampliconsFile);
        IOUtil.assertFileIsReadable(referenceSequenceFile);
        IOUtil.assertFileIsWritable(pileupCountsFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(bamFile);

        ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(
                referenceSequenceFile);

        List<Interval> amplicons = IntervalUtils.readIntervalFile(ampliconsFile);
        IntervalList intervalList = new IntervalList(referenceSequenceFileWalker.getSequenceDictionary());
        intervalList.addall(amplicons);

        SamLocusIterator locusIterator = new SamLocusIterator(reader, intervalList);

        // exclude reads that are marked as failing platform/vendor quality checks
        locusIterator.setIncludeNonPfReads(false);

        // exclude secondary alignments and reads marked as duplicates
        locusIterator.setSamFilters(Arrays.asList(new SecondaryAlignmentFilter(), new DuplicateReadFilter()));

        SamLocusAndReferenceIterator locusAndReferenceIterator = new SamLocusAndReferenceIterator(
                referenceSequenceFileWalker, locusIterator);

        BufferedWriter writer = IOUtil.openFileForBufferedWriting(pileupCountsFile);

        try {
            writeHeader(writer);

            for (SamLocusAndReferenceIterator.SAMLocusAndReference locusAndReference : locusAndReferenceIterator) {

                List<RecordAndOffset> pileup = locusAndReference.getRecordAndOffsets();

                List<RecordAndOffset> filteredPileup = PileupUtils.filterLowQualityScores(pileup, minimumBaseQuality,
                        minimumMappingQuality);

                filteredPileup = PileupUtils.filterOverlappingFragments(filteredPileup);

                writePileupCounts(writer, locusAndReference, filteredPileup);

                LocusInfo locusInfo = locusAndReference.getLocus();
                progress.record(locusInfo.getContig(), locusInfo.getPosition());
            }

            locusAndReferenceIterator.close();
            CloserUtil.close(reader);
            writer.close();
        } catch (IOException e) {
            logger.error(e);
            System.exit(1);
        }
    }

    /**
     * Write the header for the coverage/pileup output table.
     * 
     * @param writer
     * @throws IOException
     */
    private void writeHeader(BufferedWriter writer) throws IOException {
        if (id != null) {
            writer.write("ID");
            writer.write("\t");
        }
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

    private void writePileupCounts(BufferedWriter writer,
            SamLocusAndReferenceIterator.SAMLocusAndReference locusAndReference, List<RecordAndOffset> filteredPileup)
            throws IOException {

        if (id != null) {
            writer.write(id);
            writer.write("\t");
        }

        LocusInfo locusInfo = locusAndReference.getLocus();

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
