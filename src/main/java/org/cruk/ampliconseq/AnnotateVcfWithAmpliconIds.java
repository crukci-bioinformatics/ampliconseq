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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.intervals.IntervalUtils;

import htsjdk.samtools.util.CloserUtil;
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

/**
 * Utility for annotating variants in a VCF file with the identifier for the
 * overlapping amplicon.
 *
 * It is expected that each variant will overlap with one and only one amplicon.
 *
 * @author eldrid01
 */
public class AnnotateVcfWithAmpliconIds extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    private static final String AMPLICON_ATTRIBUTE = "AMPLICON";

    private File inputVcfFile;
    private File ampliconsFile;
    private File outputVcfFile;

    public static void main(String[] args) {
        AnnotateVcfWithAmpliconIds annotateVcfWithAmpliconIds = new AnnotateVcfWithAmpliconIds();
        annotateVcfWithAmpliconIds.parseCommandLineArgs(args);
        annotateVcfWithAmpliconIds.run();
    }

    @Override
    protected String getHelpDescription() {
        return "Annotates variants with the id of the overlapping amplicon";
    }

    @Override
    protected Options createOptions() {
        Options options = super.createOptions();

        Option option = new Option("i", "input", true, "VCF file (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("l", "intervals", true,
                "Amplicon intervals in which variants were called; can be in BED or Picard-style interval format (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("o", "output", true,
                "The output VCF file with variants annotated with the amplicon id (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        return options;
    }

    @Override
    protected void extractOptionValues(CommandLine commandLine) throws ParseException {
        inputVcfFile = (File) commandLine.getParsedOptionValue("input");
        ampliconsFile = (File) commandLine.getParsedOptionValue("intervals");
        outputVcfFile = (File) commandLine.getParsedOptionValue("output");
    }

    /**
     * Main run method in which VCF records are annotated with the name of the
     * overlapping amplicon.
     */
    private void run() {
        IOUtil.assertFileIsReadable(inputVcfFile);
        IOUtil.assertFileIsReadable(ampliconsFile);
        IOUtil.assertFileIsWritable(outputVcfFile);

        List<Interval> amplicons = IntervalUtils.readIntervalFile(ampliconsFile);

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

        CloserUtil.close(reader);
        writer.close();
    }
}
