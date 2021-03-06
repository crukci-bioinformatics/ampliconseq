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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;
import org.cruk.htsjdk.intervals.IntervalUtils;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Utility for extracting reads from a BAM file that correspond to targeted
 * PCR/amplicon regions.
 *
 * @author eldrid01
 */
@Command(name = "extract-amplicon-regions", versionProvider = ExtractAmpliconRegions.class, description = "\nExtracts reads from a BAM file that correspond to targeted PCR/amplicion regions.\n", mixinStandardHelpOptions = true)
public class ExtractAmpliconRegions extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    @Option(names = "--id", required = true, description = "Identifier for this dataset; if included the coverage summary will have an additional ID column (required).")
    private String id;

    @Option(names = { "-i",
            "--input" }, required = true, description = "Input BAM file which must be in coordinate sort order and indexed (required).")
    private File bamFile;

    @Option(names = { "-l",
            "--amplicon-intervals" }, required = true, description = "Amplicon intervals for which to extract matching reads; can be in BED or Picard-style interval format (required).")
    private File ampliconIntervalsFile;

    // @Option(names = "--target-intervals", required = true, description = "Target
    // intervals, can be in BED or Picard-style interval format (required).")
    // private File targetIntervalsFile;

    @Option(names = { "-d",
            "--maximum-distance" }, description = "The maximum distance of the alignment start/end to the amplicon start/end position (default: ${DEFAULT-VALUE}).")
    private int maximumDistance = 0;

    @Option(names = "--require-both-ends-anchored", description = "Require that both ends of the amplicon need to be anchored by paired end reads.")
    private boolean requireBothEndsAnchored = false;

    @Option(names = "--unmark-duplicate-reads", description = "Remove duplicate flag, if set, from SAM records.")
    private boolean unmarkDuplicateReads = false;

    @Option(names = { "-o",
            "--output" }, required = true, description = "The output BAM file containing reads that match the amplicon coordinates (required).")
    private File ampliconBamFile;

    @Option(names = { "-c",
            "--coverage" }, description = "Output coverage file summarizing read counts for each amplicon (optional).")
    private File ampliconCoverageFile;

    @Option(names = "--validation-stringency", description = "Validation stringency applied to the BAM file (default: ${DEFAULT-VALUE}).")
    private ValidationStringency validationStringency = ValidationStringency.LENIENT;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new ExtractAmpliconRegions()).execute(args);
        System.exit(exitCode);
    }

    /**
     * Main run method for extracting SAM records from a BAM file that match
     * amplicon intervals.
     */
    @Override
    public Integer call() throws Exception {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        ProgressLogger progress = new ProgressLogger(logger, 100000);

        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(ampliconIntervalsFile);
        // IOUtil.assertFileIsReadable(targetIntervalsFile);
        IOUtil.assertFileIsWritable(ampliconBamFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(validationStringency).open(bamFile);
        if (!reader.hasIndex()) {
            logger.error("No index found for input BAM file");
            return 1;
        }

        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeSAMOrBAMWriter(reader.getFileHeader(), true, ampliconBamFile);

        List<Interval> amplicons = IntervalUtils.readIntervalFile(ampliconIntervalsFile);
        // List<Interval> targets = IntervalUtils.readIntervalFile(targetIntervalsFile);

        BufferedWriter coverageWriter = null;
        if (ampliconCoverageFile != null) {
            coverageWriter = IOUtil.openFileForBufferedWriting(ampliconCoverageFile);
        }

        if (coverageWriter != null) {
            writeCoverageHeader(coverageWriter);
        }

        Map<String, Integer> ampliconReadFlags = new HashMap<>();

        // int i = 0;
        for (Interval amplicon : amplicons) {
            // Interval target = targets.get(i++);

            ampliconReadFlags.clear();

            logger.info("Extracting records for " + amplicon.toString());

            SAMRecordIterator iterator = reader.queryOverlapping(amplicon.getContig(), amplicon.getStart(),
                    amplicon.getEnd());

            // first pass over overlapping records to record which ends are
            // within an acceptable distance of the amplicon start/end
            // coordinate
            while (iterator.hasNext()) {
                SAMRecord record = iterator.next();
                if (isAcceptable(record)) {
                    addAmpliconReadFlags(record, amplicon, ampliconReadFlags);
                }
            }

            iterator.close();

            int recordCount = 0;
            int baseCount = 0;
            // int targetBaseCount = 0;

            // second pass over overlapping records in which reads or read
            // pairs consistent with this amplicon are written
            iterator = reader.queryOverlapping(amplicon.getContig(), amplicon.getStart(), amplicon.getEnd());
            while (iterator.hasNext()) {
                SAMRecord record = iterator.next();

                if (isAcceptable(record) && isAmpliconRead(record, ampliconReadFlags)) {
                    // the record is from a read or read pair that has one or
                    // both ends matching the alignment region start/end
                    // depending on whether both ends need to be anchored

                    recordCount++;

                    if (unmarkDuplicateReads) {
                        record.setDuplicateReadFlag(false);
                    }

                    writer.addAlignment(record);

                    if (coverageWriter != null) {
                        baseCount += countBasesCovered(record, amplicon);
                        // targetBaseCount += countBasesCovered(record, target);
                    }
                }

                progress.record(record);
            }
            iterator.close();

            logger.info(recordCount + " records written for " + amplicon.toString());

            if (coverageWriter != null) {
                writeCoverage(coverageWriter, amplicon, ampliconReadFlags, baseCount);
                // writeCoverage(coverageWriter, amplicon, ampliconReadFlags, baseCount, target,
                // targetBaseCount);
            }
        }

        logger.info("Writing " + ampliconBamFile.getName());
        writer.close();

        if (coverageWriter != null) {
            coverageWriter.close();
        }

        reader.close();

        logger.info("Finished");
        return 0;
    }

    /**
     * Determines whether the alignment record is acceptable for including in the
     * extracted subset matching amplicon regions.
     *
     * @param record
     * @return
     */
    private boolean isAcceptable(SAMRecord record) {
        return !(record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryAlignment()
                || record.getReadUnmappedFlag());
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
        writer.write("ID");
        writer.write("\t");
        writer.write("Amplicon");
        writer.write("\t");
        writer.write("Chromosome");
        writer.write("\t");
        writer.write("Start");
        writer.write("\t");
        writer.write("End");
        writer.write("\t");
        writer.write("Length");
        writer.write("\t");
        writer.write("Reads");
        writer.write("\t");
        writer.write("Read pairs");
        writer.write("\t");
        writer.write("Bases");
        writer.write("\t");
        writer.write("Mean coverage");
        // writer.write("\t");
        // writer.write("Target length");
        // writer.write("\t");
        // writer.write("Target bases");
        writer.write("\n");
    }

    /**
     * Write coverage metrics for a given amplicon.
     *
     * @param writer
     * @param amplicon
     * @param ampliconReadFlags
     * @param baseCount
     * @throws IOException
     */
    private void writeCoverage(BufferedWriter writer, Interval amplicon, Map<String, Integer> ampliconReadFlags,
            int baseCount) throws IOException {
        // int baseCount, Interval target, int targetBaseCount) throws IOException {

        writer.write(id);

        writer.write("\t");
        writer.write(amplicon.getName());

        writer.write("\t");
        writer.write(amplicon.getContig());
        writer.write("\t");
        writer.write(Integer.toString(amplicon.getStart()));
        writer.write("\t");
        writer.write(Integer.toString(amplicon.getEnd()));

        int length = amplicon.getEnd() - amplicon.getStart() + 1;

        writer.write("\t");
        writer.write(Integer.toString(length));

        // count reads matching at least one end of the amplicon and pairs
        // anchored to both ends
        int readCount = 0;
        int anchoredBothEndsCount = 0;
        for (int flags : ampliconReadFlags.values()) {
            if ((flags & 1) == 1) {
                readCount++;
            }
            if ((flags & 2) == 2) {
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

        writer.write("\t");
        writer.write(Integer.toString(baseCount));

        double meanCoverage = baseCount / (double) length;

        writer.write("\t");
        writer.write(Float.toString((float) meanCoverage));

        // writer.write("\t");
        // writer.write(Integer.toString(target.getEnd() - target.getStart() + 1));

        // writer.write("\t");
        // writer.write(Integer.toString(targetBaseCount));

        writer.write("\n");
    }
}
