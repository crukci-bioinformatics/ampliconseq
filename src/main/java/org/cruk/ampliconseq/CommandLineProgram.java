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

import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Abstract base class for command line programs with functions for handling
 * command-line arguments and printing help/usage on available options.
 */
public abstract class CommandLineProgram {

    /**
     * Configure command line options.
     *
     * @return an Options object
     */
    protected abstract Options createOptions();

    /**
     * Extract option values.
     *
     * @param commandLine the parsed command line object
     * @throws ParseException
     */
    protected abstract void extractOptionValues(CommandLine commandLine) throws ParseException;

    /**
     * Parse command line arguments.
     *
     * @param args the command line arguments
     */
    protected void parseCommandLineArgs(String[] args) {

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
}
