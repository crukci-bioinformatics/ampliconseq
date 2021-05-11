/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.ampliconseq;

import java.io.File;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
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
 * Utility for adding various annotations to variants in a VCF file.
 *
 * @author eldrid01
 */
@Command(name = "add-assorted-annotations-to-vcf", versionProvider = AnnotateVcfWithAmpliconIds.class, description = "\nAdds various annotations to variants in a VCF file.\n", mixinStandardHelpOptions = true)
public class AddAssortedAnnotationsToVcf extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger();

    // annotation names
    private static final String FIVE_PRIME_SEQUENCE_CONTEXT = "FivePrimeContext";
    private static final String THREE_PRIME_SEQUENCE_CONTEXT = "ThreePrimeContext";
    private static final String INDEL_LENGTH = "IndelLength";

    @Option(names = { "-i", "--input" }, required = true, description = "Input VCF file (required).")
    private File inputVcfFile;

    @Option(names = { "-r",
            "--reference-sequence" }, required = true, description = "Reference sequence FASTA file which must be indexed and have an accompanying dictionary (required).")
    private File referenceSequenceFile;

    @Option(names = { "-o",
            "--output" }, required = true, description = "Output VCF file with annotated varaints (required).")
    private File outputVcfFile;

    @Option(names = {
            "--sequence-context-length" }, description = "The number of bases of sequence context to record for both 5' and 3' context (default: ${DEFAULT-VALUE}).")
    private int sequenceContextLength = 5;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new AddAssortedAnnotationsToVcf()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public Integer call() throws Exception {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        IOUtil.assertFileIsReadable(inputVcfFile);
        IOUtil.assertFileIsReadable(referenceSequenceFile);
        IOUtil.assertFileIsWritable(outputVcfFile);

        // reference sequence
        ReferenceSequenceFile referenceSequence = ReferenceSequenceFileFactory
                .getReferenceSequenceFile(referenceSequenceFile);
        if (!referenceSequence.isIndexed()) {
            logger.error("Reference sequence FASTA file is not indexed");
            return 1;
        }

        VCFFileReader reader = new VCFFileReader(inputVcfFile, false);
        VCFHeader header = reader.getFileHeader();

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcfFile)
                .setOutputFileType(OutputType.VCF).setReferenceDictionary(header.getSequenceDictionary())
                .clearOptions();
        VariantContextWriter writer = builder.build();

        addInfoHeaderLines(header);
        writer.writeHeader(header);

        for (VariantContext variant : reader) {
            addSequenceContext(variant, referenceSequence);
            addIndelLength(variant);
            writer.add(variant);
        }

        writer.close();
        reader.close();
        referenceSequence.close();

        return 0;
    }

    /**
     * Add header lines for the added INFO fields.
     *
     * @param header the VCF header
     */
    private void addInfoHeaderLines(VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(FIVE_PRIME_SEQUENCE_CONTEXT, 1, VCFHeaderLineType.String,
                "The 5' prime sequence context."));
        header.addMetaDataLine(new VCFInfoHeaderLine(THREE_PRIME_SEQUENCE_CONTEXT, 1, VCFHeaderLineType.String,
                "The 3' prime sequence context."));
        header.addMetaDataLine(
                new VCFInfoHeaderLine(INDEL_LENGTH, 1, VCFHeaderLineType.Integer, "The length of the indel."));
    }

    /**
     * Adds INFO entries for 5' and 3' sequence context.
     *
     * @param variant the variant
     */
    private void addSequenceContext(VariantContext variant, ReferenceSequenceFile referenceSequence) {
        variant.getCommonInfo().putAttribute(FIVE_PRIME_SEQUENCE_CONTEXT,
                referenceSequence.getSubsequenceAt(variant.getContig(), variant.getStart() - sequenceContextLength,
                        variant.getStart() - 1).getBaseString());
        variant.getCommonInfo().putAttribute(THREE_PRIME_SEQUENCE_CONTEXT, referenceSequence
                .getSubsequenceAt(variant.getContig(), variant.getEnd() + 1, variant.getEnd() + sequenceContextLength)
                .getBaseString());
    }

    /**
     * Adds an INFO entry for the indel length.
     *
     * @param variant the variant
     */
    private void addIndelLength(VariantContext variant) {
        if (variant.isIndel()) {
            int length = variant.getAlternateAllele(0).length() - variant.getReference().length();
            variant.getCommonInfo().putAttribute(INDEL_LENGTH, length);
        }
    }
}
