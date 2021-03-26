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

import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.AbstractProgressLogger;

/**
 * Concrete implementation of AbstractProgressLogger using log4j for logging.
 *
 * Based on ProgressLogger from htsjdk but using log4j instead of the htsjdk Log
 * wrapper around System.err.
 *
 * @author eldrid01
 */
public class ProgressLogger extends AbstractProgressLogger {

    Logger logger = null;

    /**
     * Creates a progress logger using the given log4j logger.
     *
     * @param logger the logger
     * @param n      the frequency with which to output, e.g. every n records
     * @param verb   the verb to log, e.g. "Processed, Read, Written".
     * @param noun   the noun to use when logging, e.g. "records, variants, loci"
     */
    public ProgressLogger(Logger logger, int n, String verb, String noun) {
        super(noun, verb, n);
        this.logger = logger;
    }

    /**
     * Creates a progress logger for processing records.
     *
     * @param logger the logger
     * @param n      the frequency with which to output, e.g. every n records
     * @param noun   the noun to use when logging, e.g. "records, variants, loci"
     */
    public ProgressLogger(Logger logger, int n, String noun) {
        this(logger, n, "Processed", noun);
    }

    /**
     * Construct a progress logger for processing records.
     *
     * @param log the Log object to write outputs to
     * @param n   the frequency with which to output (i.e. every N records)
     */
    public ProgressLogger(Logger logger, int n) {
        this(logger, n, "records");
    }

    @Override
    protected void log(String... message) {
        StringBuilder sb = new StringBuilder();
        for (String s : message) {
            sb.append(s);
        }
        logger.info(sb.toString());
    }
}
