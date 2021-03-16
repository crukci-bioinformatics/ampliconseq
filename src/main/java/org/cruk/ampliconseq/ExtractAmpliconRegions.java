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
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
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
public class ExtractAmpliconRegions {
    private static final Logger logger = LogManager.getLogger();

    private File bamFile;
    private File ampliconsFile;
    private int maximumDistance;
    private boolean requireBothEndsAnchored = false;
    private boolean unmarkDuplicateReads = false;
    private File ampliconBamFile;
    private File ampliconCoverageFile;

    public static void main(String[] args) {
        ExtractAmpliconRegions extractAmpliconRegions = new ExtractAmpliconRegions();
        extractAmpliconRegions.parseCommandLineArgs(args);
        extractAmpliconRegions.run();
    }

    /**
     * Parse command line arguments.
     *
     * @param args
     */
    private void parseCommandLineArgs(String[] args) {

        Options options = createOptions();

        checkForHelp(args, options);

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(options, args);

            extractOptionValues(commandLine);

        } catch (ParseException e) {
            System.err.println("Error parsing command-line arguments");
            System.err.println();
            System.err.println(e.getMessage());
            printUsage(options, System.err);
            System.exit(1);
        }
    }

    /**
     * Configure command line options.
     *
     * @return
     */
    private Options createOptions() {
        Options options = new Options();

        Option option = new Option("h", "help", false, "Print command line options");
        options.addOption(option);

        option = new Option("i", "bam", true, "BAM input file(s), which must be in coordinate sort order and indexed.");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("a", "amplicons", true,
                "The file containing the set of amplicons from which to create non-overlapping sets; can be in BED or Picard-style interval format.");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option(null, "maximum-distance", true,
                "The maximum distance of the alignment start/end to the amplicon start/end position.");
        option.setRequired(false);
        option.setType(Number.class);
        options.addOption(option);

        option = new Option(null, "require-both-ends-anchored", false,
                "Set if both ends of the amplicon need to be anchored by paired end reads.");
        option.setRequired(false);
        options.addOption(option);

        option = new Option(null, "unmark-duplicate-reads", false, "Remove duplicate flag, if set, from reads.");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("o", "amplicon-bam", true,
                "The output BAM file containing reads for the given set of amplicons.");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("c", "coverage", true,
                "Optional output coverage file summarizing read counts for each amplicon.");
        option.setRequired(false);
        option.setType(File.class);
        options.addOption(option);

        return options;
    }

    /**
     * Extract option values.
     *
     * @param commandLine
     * @throws ParseException
     */
    private void extractOptionValues(CommandLine commandLine) throws ParseException {
        bamFile = (File) commandLine.getParsedOptionValue("bam");
        ampliconsFile = (File) commandLine.getParsedOptionValue("amplicons");
        if (commandLine.hasOption("maximum-distance")) {
            maximumDistance = ((Number) commandLine.getParsedOptionValue("maximum-distance")).intValue();
        }
        requireBothEndsAnchored = commandLine.hasOption("require-both-ends-anchored");
        unmarkDuplicateReads = commandLine.hasOption("unmark-duplicate-reads");
        ampliconBamFile = (File) commandLine.getParsedOptionValue("amplicon-bam");
        ampliconCoverageFile = (File) commandLine.getParsedOptionValue("coverage");
    }

    /**
     * Workaround to allow for a help message to be displayed to stdout without
     * error or warning messages about missing arguments.
     * 
     * @param args
     * @param options
     */
    private void checkForHelp(String[] args, Options options) {
        if (!options.hasLongOption("help"))
            return;

        // create a copy of the options with none set to be required
        Options checkForHelpOptions = new Options();
        for (Option option : options.getOptions()) {
            checkForHelpOptions.addOption(option.getOpt(), option.getLongOpt(), option.hasArg(),
                    option.getDescription());
        }

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(checkForHelpOptions, args, true);
            if (commandLine.hasOption("help")) {
                printUsage(options, System.out);
                System.exit(0);
            }
        } catch (ParseException e) {
        }
    }

    /**
     * Print help/usage.
     *
     * @param options
     * @param stream
     */
    private void printUsage(Options options, PrintStream stream) {
        PrintWriter out = new PrintWriter(new OutputStreamWriter(stream));
        out.println();
        HelpFormatter helpFormatter = new HelpFormatter();
        String syntax = "java " + getClass().getName() + " [options]";
        String description = "\nExtracts reads from a BAM file that correspond to targeted PCR/amplicion regions.\n\n";
        helpFormatter.printHelp(out, 80, syntax, description, options, 4, 8, "", false);
        out.println();
        out.flush();
    }

    private void run() {
        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(ampliconsFile);
        IOUtil.assertFileIsWritable(ampliconBamFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(bamFile);
        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeSAMOrBAMWriter(reader.getFileHeader(), false, ampliconBamFile);

        BufferedWriter coverageWriter = null;
        if (ampliconCoverageFile != null)
            coverageWriter = IOUtil.openFileForBufferedWriting(ampliconCoverageFile);

        try {
            List<Interval> amplicons = IntervalUtils.readIntervalFile(ampliconsFile);

            writeCoverageHeader(coverageWriter);

            Map<String, Integer> observedLookup = new HashMap<String, Integer>();

            for (Interval amplicon : amplicons) {
                observedLookup.clear();

                logger.info("Extracting records for " + amplicon.toString());

                SAMRecordIterator iterator = reader.queryOverlapping(amplicon.getContig(), amplicon.getStart(),
                        amplicon.getEnd());

                // first pass over overlapping records to record which ends are
                // within an acceptable distance of the amplicon start/end
                // coordinate
                while (iterator.hasNext()) {
                    SAMRecord record = iterator.next();

                    if (record.isSecondaryAlignment() || record.getReadUnmappedFlag())
                        continue;

                    String readName = record.getReadName();
                    int observedFlag = observedLookup.getOrDefault(readName, 0);

                    int alignmentStartDistance = Math.abs(record.getAlignmentStart() - amplicon.getStart());
                    int unclippedStartDistance = Math.abs(record.getUnclippedStart() - amplicon.getStart());
                    int alignmentEndDistance = Math.abs(record.getAlignmentEnd() - amplicon.getEnd());
                    int unclippedEndDistance = Math.abs(record.getUnclippedEnd() - amplicon.getEnd());

                    if (alignmentStartDistance <= maximumDistance && unclippedStartDistance <= maximumDistance
                            && !record.getReadNegativeStrandFlag()) {
                        observedFlag |= 1;
                        if (record.getReadPairedFlag())
                            observedFlag |= (record.getFirstOfPairFlag() ? 4 : 8);
                    } else if (alignmentEndDistance <= maximumDistance && unclippedEndDistance <= maximumDistance
                            && record.getReadNegativeStrandFlag()) {
                        observedFlag |= 2;
                        if (record.getReadPairedFlag())
                            observedFlag |= (record.getFirstOfPairFlag() ? 4 : 8);
                    }

                    if (observedFlag != 0)
                        observedLookup.put(readName, observedFlag);
                }

                iterator.close();

                int recordCount = 0;
                int baseCount = 0;

                // second pass over overlapping records in which reads or read
                // pairs consistent with this amplicon are written
                iterator = reader.queryOverlapping(amplicon.getContig(), amplicon.getStart(), amplicon.getEnd());
                while (iterator.hasNext()) {
                    SAMRecord record = iterator.next();
                    if (record.isSecondaryAlignment())
                        continue;
                    String readName = record.getReadName();
                    if (observedLookup.containsKey(readName)) {
                        int observedFlag = observedLookup.get(readName);
                        if (!requireBothEndsAnchored || observedFlag == 15) {
                            recordCount++;

                            if (unmarkDuplicateReads)
                                record.setDuplicateReadFlag(false);

                            writer.addAlignment(record);

                            baseCount += countBasesCovered(record, amplicon.getStart(), amplicon.getEnd());
                        }
                    }
                }
                iterator.close();

                logger.info(recordCount + " records written");

                if (coverageWriter != null) {
                    // pair count where both ends anchored to different ends of the amplicon
                    int anchoredBothEndsCount = 0;
                    for (int observedFlag : observedLookup.values()) {
                        if (observedFlag == 15)
                            anchoredBothEndsCount++;
                    }
                    writeCoverage(coverageWriter, amplicon, baseCount, observedLookup.size(), anchoredBothEndsCount);
                }
            }

            CloserUtil.close(reader);
            writer.close();
            coverageWriter.close();

        } catch (IOException e) {
            logger.error(e);
            System.exit(1);
        }
    }

    /**
     * Count bases covered within the given genomic interval.
     *
     * @param record
     * @param start
     * @param end
     * @return
     */
    private int countBasesCovered(SAMRecord record, int start, int end) {
        int count = 0;
        if (!record.getReadFailsVendorQualityCheckFlag() && !record.getReadUnmappedFlag()
                && record.getMappingQuality() != 0 && !record.getDuplicateReadFlag()) {
            for (AlignmentBlock block : record.getAlignmentBlocks()) {
                int last = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                for (int pos = block.getReferenceStart(); pos <= last; ++pos) {
                    if (pos >= start && pos <= end) {
                        count++;
                    }
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
        if (writer != null) {
            writer.write("CHROM\tSTART\tEND\tLENGTH\tNAME\tMEAN_COVERAGE\tREADS\tREAD_PAIRS\n");
        }
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
    private void writeCoverage(BufferedWriter writer, Interval amplicon, int baseCount, int pairCount,
            int anchoredBothEndsCount) throws IOException {
        if (writer != null) {

            int length = amplicon.getEnd() - amplicon.getStart() + 1;
            double meanCoverage = baseCount / (double) (amplicon.getEnd() - amplicon.getStart() + 1);

            writer.write(amplicon.getContig());
            writer.write("\t");
            writer.write(Integer.toString(amplicon.getStart()));
            writer.write("\t");
            writer.write(Integer.toString(amplicon.getEnd()));
            writer.write("\t");
            writer.write(Integer.toString(length));
            writer.write("\t");
            writer.write(amplicon.getName());
            writer.write("\t");
            writer.write(Float.toString((float) meanCoverage));
            writer.write("\t");
            writer.write(Integer.toString(pairCount));
            writer.write("\t");
            writer.write(Integer.toString(anchoredBothEndsCount));
            writer.write("\n");
        }
    }
}
