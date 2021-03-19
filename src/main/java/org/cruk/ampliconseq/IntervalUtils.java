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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

/**
 * Utility methods for reading intervals from BED and Picard-style files.
 *
 * @author eldrid01
 */
public class IntervalUtils {

    /**
     * Reads intervals from a BED file or Picard-style interval list format
     * depending on file extension.
     *
     * @param intervalFile the BED or interval list file.
     * @return a list of intervals
     */
    public static List<Interval> readIntervalFile(File intervalFile) {
        if (intervalFile.getName().toLowerCase().endsWith(".bed")) {
            return IntervalUtils.readBEDFile(intervalFile);
        } else {
            return IntervalUtils.readIntervalListFile(intervalFile);
        }
    }

    /**
     * Reads a Picard-style interval list file (5 column tabular format with
     * sequence dictionary header).
     *
     * Note that an unmodifiable list of Interval objects is returned by the
     * IntervalList created by IntervalList.fromFile(). Here the contents are copied
     * to a regular ArrayList.
     *
     * @param intervalListFile the Picard-style interval list file.
     * @return a list of intervals
     */
    public static List<Interval> readIntervalListFile(File intervalListFile) {
        IntervalList intervalList = IntervalList.fromFile(intervalListFile);
        // note that an unmodifiable list is returned by intervalList.getIntervals()
        List<Interval> intervals = new ArrayList<Interval>();
        intervals.addAll(intervalList.getIntervals());
        return intervals;

    }

    /**
     * Reads intervals from a BED file.
     *
     * @param bedFile the BED format interval file.
     * @return a list of intervals
     */
    public static List<Interval> readBEDFile(File bedFile) {
        FeatureReader<BEDFeature> reader = AbstractFeatureReader.getFeatureReader(bedFile.getAbsolutePath(),
                new BEDCodec(), false);
        try {
            CloseableTribbleIterator<BEDFeature> iterator = reader.iterator();
            List<Interval> intervals = new ArrayList<Interval>();
            while (iterator.hasNext()) {
                BEDFeature feature = iterator.next();
                Interval interval = new Interval(feature.getContig(), feature.getStart(), feature.getEnd(),
                        feature.getStrand() == Strand.NEGATIVE, feature.getName());
                intervals.add(interval);
            }
            return intervals;
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    throw new RuntimeIOException(e);
                }
            }
        }
    }
}
