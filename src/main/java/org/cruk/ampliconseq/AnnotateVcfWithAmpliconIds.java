/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.ampliconseq;

import java.io.File;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.intervals.IntervalUtils;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Utility for annotating variants in a VCF file with the identifier for the
 * overlapping amplicon.
 *
 * It is expected that each variant will overlap with one and only one amplicon.
 *
 * @author eldrid01
 */
@Command(name = "annotate-vcf-with-amplicon-ids", versionProvider = AnnotateVcfWithAmpliconIds.class, description = "\nAnnotates variants with the id of the overlapping amplicon.\n", mixinStandardHelpOptions = true)
public class AnnotateVcfWithAmpliconIds extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    private static final String AMPLICON_ATTRIBUTE = "AMPLICON";

    @Option(names = { "-i", "--input" }, required = true, description = "Input VCF file (required).")
    private File inputVcfFile;

    @Option(names = { "-l",
            "--amplicon-intervals" }, required = true, description = "Amplicon intervals in which variants were called; can be in BED or Picard-style interval format (required).")
    private File ampliconIntervalsFile;

    @Option(names = { "-o",
            "--output" }, required = true, description = "Output VCF file with variants annotated with the amplicon id (required).")
    private File outputVcfFile;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new AnnotateVcfWithAmpliconIds()).execute(args);
        System.exit(exitCode);
    }

    /**
     * Main run method in which VCF records are annotated with the name of the
     * overlapping amplicon.
     */
    @Override
    public Integer call() throws Exception {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        IOUtil.assertFileIsReadable(inputVcfFile);
        IOUtil.assertFileIsReadable(ampliconIntervalsFile);
        IOUtil.assertFileIsWritable(outputVcfFile);

        List<Interval> amplicons = IntervalUtils.readIntervalFile(ampliconIntervalsFile);

        VCFFileReader reader = new VCFFileReader(inputVcfFile, false);
        VCFHeader header = reader.getFileHeader();

        header.addMetaDataLine(
                new VCFInfoHeaderLine(AMPLICON_ATTRIBUTE, 1, VCFHeaderLineType.String, "The amplicon identifier."));

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcfFile)
                .setOutputFileType(OutputType.VCF).setReferenceDictionary(header.getSequenceDictionary())
                .clearOptions();
        VariantContextWriter writer = builder.build();

        writer.writeHeader(header);

        for (VariantContext variant : reader) {
            String ampliconId = null;

            for (Interval amplicon : amplicons) {
                if (variant.getContig().equals(amplicon.getContig()) && variant.getStart() <= amplicon.getEnd()
                        && variant.getEnd() >= amplicon.getStart()) {
                    if (ampliconId == null) {
                        ampliconId = amplicon.getName();
                    } else {
                        logger.warn("Multiple amplicons found for variant " + variant);
                    }
                }
            }

            if (ampliconId == null) {
                // VarDict has been observed to call variants just outside a target region
                logger.warn("No overlapping amplicon for variant " + variant);
            } else {
                variant.getCommonInfo().putAttribute(AMPLICON_ATTRIBUTE, ampliconId);
            }

            writer.add(variant);
        }

        writer.close();

        reader.close();

        logger.info("Finished");
        return 0;
    }

    // annotation names
    private static final String INDEL_LENGTH = "IndelLength";

    /**
     * Add header lines for the added INFO fields.
     *
     * @param header the VCF header
     */
    private void addInfoHeaderLines(VCFHeader header) {
        header.addMetaDataLine(
                new VCFInfoHeaderLine(INDEL_LENGTH, 1, VCFHeaderLineType.Integer, "The length of the indel."));
    }

    private void addIndelLength(VariantContext variant) {
        if (variant.isIndel()) {
            
        }
    }
}
