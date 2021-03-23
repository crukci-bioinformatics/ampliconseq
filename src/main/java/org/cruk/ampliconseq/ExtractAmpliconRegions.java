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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;

/**
 * Utility for extracting reads from a BAM file that correspond to targeted
 * PCR/amplicon regions.
 *
 * @author eldrid01
 */
public class ExtractAmpliconRegions extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    private String id;
    private File bamFile;
    private File ampliconsFile;
    private int maximumDistance = 0;
    private boolean requireBothEndsAnchored = false;
    private boolean unmarkDuplicateReads = false;
    private File ampliconBamFile;
    private File ampliconCoverageFile;

    public static void main(String[] args) {
        ExtractAmpliconRegions extractAmpliconRegions = new ExtractAmpliconRegions();
        extractAmpliconRegions.parseCommandLineArgs(args);
        extractAmpliconRegions.run();
    }

    @Override
    protected String getHelpDescription() {
        return "Extracts reads from a BAM file that correspond to targeted PCR/amplicion regions";
    }

    @Override
    protected Options createOptions() {
        Options options = new Options();

        Option option = new Option("h", "help", false, "Print command line options");
        options.addOption(option);

        option = new Option(null, "id", true,
                "Identifier for this dataset; if included the coverage summary will have an additional ID column (optional");
        options.addOption(option);

        option = new Option("i", "bam", true,
                "BAM input file which must be in coordinate sort order and indexed (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("l", "intervals", true,
                "Amplicon intervals for which to extract matching reads; can be in BED or Picard-style interval format (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("d", "maximum-distance", true,
                "The maximum distance of the alignment start/end to the amplicon start/end position (default: "
                        + maximumDistance + ")");
        option.setType(Number.class);
        options.addOption(option);

        option = new Option(null, "require-both-ends-anchored", false,
                "Set if both ends of the amplicon need to be anchored by paired end reads");
        options.addOption(option);

        option = new Option(null, "unmark-duplicate-reads", false, "Remove duplicate flag, if set, from reads");
        options.addOption(option);

        option = new Option("o", "amplicon-bam", true,
                "The output BAM file containing reads that match the amplicon coordinates (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("c", "coverage", true,
                "Output coverage file summarizing read counts for each amplicon (optional)");
        option.setType(File.class);
        options.addOption(option);

        return options;
    }

    @Override
    protected void extractOptionValues(CommandLine commandLine) throws ParseException {
        id = commandLine.getOptionValue("id");
        bamFile = (File) commandLine.getParsedOptionValue("bam");
        ampliconsFile = (File) commandLine.getParsedOptionValue("intervals");
        if (commandLine.hasOption("maximum-distance")) {
            maximumDistance = ((Number) commandLine.getParsedOptionValue("maximum-distance")).intValue();
        }
        requireBothEndsAnchored = commandLine.hasOption("require-both-ends-anchored");
        unmarkDuplicateReads = commandLine.hasOption("unmark-duplicate-reads");
        ampliconBamFile = (File) commandLine.getParsedOptionValue("amplicon-bam");
        ampliconCoverageFile = (File) commandLine.getParsedOptionValue("coverage");
    }

    /**
     * Main run method for extracting SAM records from a BAM file that match
     * amplicon intervals.
     */
    private void run() {

        ProgressLogger progress = new ProgressLogger(logger, 100000);

        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(ampliconsFile);
        IOUtil.assertFileIsWritable(ampliconBamFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(bamFile);
        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeSAMOrBAMWriter(reader.getFileHeader(), false, ampliconBamFile);

        List<Interval> amplicons = IntervalUtils.readIntervalFile(ampliconsFile);

        BufferedWriter coverageWriter = null;
        if (ampliconCoverageFile != null) {
            coverageWriter = IOUtil.openFileForBufferedWriting(ampliconCoverageFile);
        }

        try {
            if (coverageWriter != null) {
                writeCoverageHeader(coverageWriter);
            }

            Map<String, Integer> ampliconReadFlags = new HashMap<String, Integer>();

            for (Interval amplicon : amplicons) {
                ampliconReadFlags.clear();

                logger.info("Extracting records for " + amplicon.toString());

                SAMRecordIterator iterator = reader.queryOverlapping(amplicon.getContig(), amplicon.getStart(),
                        amplicon.getEnd());

                // first pass over overlapping records to record which ends are
                // within an acceptable distance of the amplicon start/end
                // coordinate
                while (iterator.hasNext()) {
                    SAMRecord record = iterator.next();
                    addAmpliconReadFlags(record, amplicon, ampliconReadFlags);
                }

                iterator.close();

                int recordCount = 0;
                int baseCount = 0;

                // second pass over overlapping records in which reads or read
                // pairs consistent with this amplicon are written
                iterator = reader.queryOverlapping(amplicon.getContig(), amplicon.getStart(), amplicon.getEnd());
                while (iterator.hasNext()) {
                    SAMRecord record = iterator.next();

                    if (isAmpliconRead(record, ampliconReadFlags)) {
                        recordCount++;

                        if (unmarkDuplicateReads) {
                            record.setDuplicateReadFlag(false);
                        }

                        writer.addAlignment(record);

                        if (coverageWriter != null) {
                            baseCount += countBasesCovered(record, amplicon);
                        }
                    }

                    progress.record(record);
                }
                iterator.close();

                logger.info(recordCount + " records written for " + amplicon.toString());

                if (coverageWriter != null) {
                    writeCoverage(coverageWriter, amplicon, ampliconReadFlags, baseCount);
                }
            }

            logger.info("Writing " + ampliconBamFile.getName());
            writer.close();

            CloserUtil.close(reader);

            if (coverageWriter != null) {
                coverageWriter.close();
            }

        } catch (IOException e) {
            logger.error(e);
            System.exit(1);
        }
    }

    /**
     * Sets flags for a read that overlaps a given amplicon.
     * 
     * @param record            the SAM record
     * @param amplicon          the amplicon interval
     * @param ampliconReadFlags a lookup of bit flags for reads overlapping the
     *                          amplicon
     */
    private void addAmpliconReadFlags(SAMRecord record, Interval amplicon, Map<String, Integer> ampliconReadFlags) {
        if (record.isSecondaryAlignment() || record.getReadUnmappedFlag()) {
            return;
        }

        String name = record.getReadName();
        int flags = ampliconReadFlags.getOrDefault(name, 0);

        int matchFlag = 0;
        if (record.getReadNegativeStrandFlag()) {
            if (Math.abs(record.getAlignmentEnd() - amplicon.getEnd()) <= maximumDistance
                    && Math.abs(record.getUnclippedEnd() - amplicon.getEnd()) <= maximumDistance) {
                matchFlag = 2;
            }
        } else {
            if (Math.abs(record.getAlignmentStart() - amplicon.getStart()) <= maximumDistance
                    && Math.abs(record.getUnclippedStart() - amplicon.getStart()) <= maximumDistance) {
                matchFlag = 1;
            }
        }

        if (matchFlag == 0) {
            return;
        }

        int readInPairFlag = 0;
        if (record.getReadPairedFlag()) {
            readInPairFlag = record.getFirstOfPairFlag() ? 4 : 8;
        }

        flags |= matchFlag;
        flags |= readInPairFlag;

        ampliconReadFlags.put(name, flags);
    }

    /**
     * Checks the bit flags to see if the given SAM record comes from a read matched
     * to the current amplicon.
     *
     * @param record            the SAM record
     * @param ampliconReadFlags a lookup of bit flags for reads overlapping the
     *                          amplicon
     * @return <code>true</code> if the record is for a read that matches the
     *         current amplicon
     */
    private boolean isAmpliconRead(SAMRecord record, Map<String, Integer> ampliconReadFlags) {
        if (record.isSecondaryAlignment() || record.getReadUnmappedFlag()) {
            return false;
        }

        String name = record.getReadName();
        int flags = ampliconReadFlags.getOrDefault(name, 0);

        if (requireBothEndsAnchored) {
            return flags == 15;
        } else {
            return (flags & 3) != 0;
        }
    }

    /**
     * Count bases covered within the given genomic interval.
     *
     * @param record   the SAM record
     * @param amplicon the amplicon interval
     * @return
     */
    private int countBasesCovered(SAMRecord record, Interval amplicon) {
        if (record.getReadFailsVendorQualityCheckFlag() || record.getReadUnmappedFlag() || record.getDuplicateReadFlag()
                || record.getMappingQuality() == 0)
            return 0;

        int count = 0;
        for (AlignmentBlock block : record.getAlignmentBlocks()) {
            int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
            for (int pos = block.getReferenceStart(); pos <= end; ++pos) {
                if (pos >= amplicon.getStart() && pos <= amplicon.getEnd()) {
                    count++;
                }
            }
        }
        return count;
    }

    /**
     * Write the coverage metrics header.
     * 
     * @param writer
     * @throws IOException
     */
    private void writeCoverageHeader(BufferedWriter writer) throws IOException {
        if (id != null) {
            writer.write("ID");
            writer.write("\t");
        }
        writer.write("Chromosome");
        writer.write("\t");
        writer.write("Start");
        writer.write("\t");
        writer.write("End");
        writer.write("\t");
        writer.write("Length");
        writer.write("\t");
        writer.write("Name");
        writer.write("\t");
        writer.write("Mean coverage");
        writer.write("\t");
        writer.write("Reads");
        writer.write("\t");
        writer.write("Read pairs");
        writer.write("\n");
    }

    /**
     * Write coverage metrics for a given amplicon.
     * 
     * @param writer
     * @param amplicon
     * @param baseCount
     * @param pairCount
     * @param anchoredBothEndsCount
     * @throws IOException
     */
    private void writeCoverage(BufferedWriter writer, Interval amplicon, Map<String, Integer> ampliconReadFlags,
            int baseCount) throws IOException {

        if (id != null) {
            writer.write(id);
            writer.write("\t");
        }

        writer.write(amplicon.getContig());
        writer.write("\t");
        writer.write(Integer.toString(amplicon.getStart()));
        writer.write("\t");
        writer.write(Integer.toString(amplicon.getEnd()));

        int length = amplicon.getEnd() - amplicon.getStart() + 1;

        writer.write("\t");
        writer.write(Integer.toString(length));

        writer.write("\t");
        writer.write(amplicon.getName());

        double meanCoverage = baseCount / (double) length;

        writer.write("\t");
        writer.write(Float.toString((float) meanCoverage));

        // count reads matching at least one end of the amplicon and pairs
        // anchored to both ends
        int readCount = 0;
        int anchoredBothEndsCount = 0;
        for (int flags : ampliconReadFlags.values()) {
            if (flags != 0) {
                readCount++;
            }
            if (flags == 15) {
                anchoredBothEndsCount++;
            }
        }

        writer.write("\t");
        writer.write(Integer.toString(readCount));
        writer.write("\t");
        writer.write(Integer.toString(anchoredBothEndsCount));

        writer.write("\n");
    }
}
