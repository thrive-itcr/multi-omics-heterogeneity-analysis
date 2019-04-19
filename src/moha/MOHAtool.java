package moha;

/**
 * Copyright (c) General Electric Company, 2018.
 * All rights reserved.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.swing.table.DefaultTableModel;


/**
 * The Multi-Omics Heterogeneity Analysis (MOHA) tool.
 * <p>
 * Copyright (c) General Electric Company, 2018.<br />
 * All rights reserved.
 * 
 * @version 1.0.0 
 * 
 * @author GE Global Research
 * <p>
 * Repository Name: moha-heterogeneity-metrics<br />
 * Repository Description: Compute heterogeneity metrics using Multi-Omics Heterogeneity Analysis (MOHA) tool
 * <p>
 * The MOHA algorithm and approach is described in the publication:<br />
 * Graf JF, Zavodszky MI (2017) Characterizing the heterogeneity of tumor tissues from spatially resolved molecular measures. PLoS ONE 12: e0188878.<br />
 * <a href="https://www.ncbi.nlm.nih.gov/pubmed/29190747">https://www.ncbi.nlm.nih.gov/pubmed/29190747</a>
 * <p>
 * The MOHA tool runs from the command line and performs four major functions.
 * <ol>
 * <li>Computes thresholds for cell measures and saves them to a file.
 * <br />Usage:<p>
 * java -jar MOHAtool.jar  [-options] -computeThresholds="input Study or Cell Measurement file" -thresholdFile="output Threshold file"
 * <p>
 * Example using the test data.
 * <p>
 * java -jar MOHAtool.jar  -computeThresholds=test/quant_006.csv -thresholdFile=test/myThresholds.txt
 * <p>
 * </li>
 * <li>Compute cell states from inputted cell measure files, and an existing threshold file. The cell states are outputted and saved to a cell marker state file. A marker index file that relates markers to genes is also generated and saved.
 * <br />Usage:<p>
 * java -jar MOHAtool.jar  [-options] -computeCellStates="input Study or Cell Measurement file" -thresholdFile="input Threshold file" -cmsfile="output Cell Marker State file" -mifile="output Marker Index file"
 * <p>
 * Example using the test data.
 * <p>
 * java -jar MOHAtool.jar  -computeCellStates=test/quant_006.csv -thresholdFile=test/myThresholds.txt -cmsfile=test/myOutput.MarkerStates.txt -mifile=test/myOutput.MarkerIndex.txt
 * <p>
 * </li>
 * <li>Computes the heterogeneity metrics for multiple samples of cells (FOVs). The marker combinations that are used come from inputted gene sets (e.g. AKT pathway).
 * <br />Usage:<p>
 * java -jar MOHAtool.jar  [-options] -computeGeneSetHeterogeneity="input Gene Set file" -cmsfile="input Cell Marker State file" -mifile="input Marker Index file" -outputFile="output file"
 * <p>
 * Two examples using the test data.
 * <p>
 * java -jar MOHAtool.jar  -computeGeneSetHeterogeneity=test/exampleMarkerSetFile.txt -cmsfile=test/myOutput.MarkerStates.txt -mifile=test/myOutput.MarkerIndex.txt -outputFile=test/myOutput.moha.txt
 * <p>
 * java -jar MOHAtool.jar  -computeGeneSetHeterogeneity=test/crcStudy.GeneSets.txt -cmsfile=test/crcStudy.MarkerStates.txt -mifile=test/crcStudy.MarkerIndex.txt -outputFile=test/myOutput.moha.txt
 * <p>
 * </li>
 * <li>Computes the heterogeneity metrics for an inputted cell state file representing a single sample (FOV or position on a slide).
 * <br />Usage:<p>
 * java -jar MOHAtool.jar  [-options] -computeHeterogeneity="input Cell Marker State file" -outputFile="output file"
 * <p>
 * Example using the test data.
 * <p>
 * java -jar MOHAtool.jar  -computeHeterogeneity=test/AGA_260_3_AKT.MarkerStates.txt -outputFile=test/myOutput.moha.txt -append=false
 * <p>
 * </li>
 * </ol>
 * <p>
 * The MOHA tool works with multiple types of text files that are either tab or comma delimited. 
 * The first line of each of these text files must contain the data column headings. Furthermore, 
 * for each file type there are a few columns that are always required while others may be optional.
 * <ol>
 * <li>Cell measurement files are used as input to the MOHA tool.
 * These files contain both cell spatial information and biomarker measurements. Each row in the file is an individual cell.
 * The file can contain multiple biomarker measurement columns with the name of the column following a format of Biomarker_Locataion_Metric.
 * For example, if the biomarker is AKT and its Median intensity is measured for the whole Cell, then the column name would be AKT_Cell_Median
 * Two columns are required to provide the spatial X and Y coordinates for each cell.
 * <ul><li>Cell_Center_X</li><li>Cell_Center_Y</li></ul>
 * Although not required, the Cell_ID and Cell_Area are two columns that will be used if provided in a cell measurement file. 
 * <ul><li>Cell_ID</li><li>Cell_Area</li></ul>
 * <p></li>
 * <li>Threshold files contain values that are used to convert continuous measurements into ordinal marker state values for each Biomarker_Locataion_Metric.
 * The MOHA tool creates a Threshold file as output when running the computeThresholds method. A threshold file is required as input when running the computeCellStates method
 * <p></li>
 * <li>Cell marker state file is created as output when running the computeCellStates method.
 * This file type is used as input when running the computeHeterogeneity and computeGeneSetHeterogeneity methods.
 * The required columns for a cell marker state file (cmsFile) include:
 * <ul><li>Cell_Center_X</li><li>Cell_Center_Y</li><li>Marker_States</li></ul>
 * Although not required the following optional columns will be used if provided in the file. 
 * <ul><li>Sample_ID</li><li>Slide_ID</li><li>Position_ID</li><li>Cell_ID</li><li>Cell_Area</li><li>Cell_Radius</li></ul>
 * <p></li>
 * <li>Cell marker index file is created as output when running the computeCellStates method.
 * This file type is used as input when running the computeGeneSetHeterogeneity method.
 * The required columns for a cell marker index file (miFile) include:
 * <ul><li>MarkerName</li><li>MarkerTargets</li></ul>
 * Although not required the following optional columns will be used if provided in the file. 
 * <ul><li>MarkerIndx</li><li>LocationName</li></ul>
 * <p></li>
 * <li>Study file can be used to input a list of cell measurement files when running either the computeThresholds or computeCellStates method.
 * The required columns for a study file include:
 * <ul><li>Sample_ID</li><li>DATA_FILENAME</li></ul>
 * <p></li>
 * <li>Gene set file is used as input when running the computeGeneSetHeterogeneity method.
 * The required columns for a Gene Set (or Marker Set) file include:
 * <ul><li>markerSetName</li><li>markers</li></ul>
 * <p></li>
 * <li>Diversity file is created when running either the computeHeterogeneity or computeGeneSetHeterogeneity methods.
 * This file contains the computed heterogeneity metrics for each collection of cells or samples.
 * <p></li>
 * </ol>
 */



public class MOHAtool {

	protected final static String versionNumber = "1.0.0";

	protected static final String option_computeThresholds = "computeThresholds";
	protected static final String option_computeCellStates = "computeCellStates";
	protected static final String option_computeGeneSetHeterogeneity = "computeGeneSetHeterogeneity";
	protected static final String option_computeHeterogeneity = "computeHeterogeneity";

	protected static final String option_runTestComputeThresholds = "runTestComputeThresholds";
	protected static final String option_runTestComputeCellStates = "runTestComputeCellStates";
	protected static final String option_runTestGeneSetHeterogeneity = "runTestGeneSetHeterogeneity";
	protected static final String option_runTestHeterogeneity = "runTestHeterogeneity";

	protected static final String option_version = "version";
	protected static final String option_help = "help";
	protected static final String option_h = "h";
	protected static final String option_questionMark = "?";
	protected static final String option_verbose = "verbose";
	protected static final String option_outputDir = "outputDir";
	protected static final String option_outputFile = "outputFile";
	protected static final String option_append = "append";
	protected static final String option_cdf = "cdf";
	protected static final String option_nStateModel = "nStateModel";
	protected static final String option_SampleID = "SampleID";
	protected static final String option_GeneSetName = "GeneSetName";
	protected static final String option_MaxNumCellStates = "MaxNumCellStates";
	protected static final String option_MaxNumCells = "MaxNumCells";
	
	protected static final String option_cmsFile = "cmsfile";   // cell marker state file
	protected static final String option_miFile = "mifile";    // marker index file
	protected static final String option_thresholdFile = "thresholdFile";
	
	protected static final String option_colName_Sample_ID = "colName_Sample_ID";
	protected static final String option_colName_DATA_FILENAME = "colName_DATA_FILENAME";
	protected static final String option_colName_Slide_ID = "colName_Slide_ID";
	protected static final String option_colName_Position_ID = "colName_Position_ID";
	protected static final String option_colName_Cell_ID = "colName_Cell_ID";
	protected static final String option_colName_Cell_Center_X = "colName_Cell_Center_X";
	protected static final String option_colName_Cell_Center_Y = "colName_Cell_Center_Y";
	protected static final String option_colName_Cell_Area = "colName_Cell_Area";
	protected static final String option_colName_Cell_Radius = "colName_Cell_Radius";
	protected static final String option_colName_Marker_States = "colName_Marker_States";
	protected static final String option_biomarkerColNameTag = "biomarkerColNameTag";
	
	
	protected DecimalFormat df4 = new DecimalFormat("0.0000");
	
	protected boolean verbose = false;
	
	
	/**
	 * The class main method.
	 * <p>
	 * To test the methods, change working directory to location of "MOHAtool.jar" that should have a "test" subdirectory containing files that end with "_Validation.txt" and
	 * then run the following four command lines from the console.
	 * <p>
	 * java -jar MOHAtool.jar -runTestComputeThresholds=test/quant_006.csv.thresholds_Validation.txt
	 * <p>
	 * java -jar MOHAtool.jar -runTestComputeCellStates=test/quant_006.csv.MarkerStates_Validation.txt
	 * <p>
	 * java -jar MOHAtool.jar -runTestHeterogeneity=test/AGA_260_3_AKT_Validation.txt
	 * <p>
	 * java -jar MOHAtool.jar -runTestGeneSetHeterogeneity=test/crcStudy_Validation.txt
	 * <p>
	 * @param args arguments of command line after parsing by space
	 */
    public static void main(String[] args) {

	    MOHAtool mohaTool = new MOHAtool();
	    mohaTool.runTool(args);

    }
   
    
    
	/**
	 * This method is called by the class main upon starting tool execution 
	 */
	protected void runTool(String[] args) {
		
    	boolean RUN_DEVELOPMENT_MODE = false;
    	if (RUN_DEVELOPMENT_MODE) {
    		System.out.println("--------------------------------------------");
    		System.out.println("Running in development mode");
    		System.out.println("--------------------------------------------");
    		

//    		java -jar MOHAtool.jar -runTestComputeThresholds=test/quant_006.csv.thresholds_Validation.txt
//    		args = new String[1];
//    		args[0] = "-" + option_runTestComputeThresholds + "=test" + File.separator + "quant_006.csv.thresholds_Validation.txt";


//    		java -jar MOHAtool.jar -runTestComputeCellStates=test/quant_006.csv.MarkerStates_Validation.txt
//    		args = new String[1];
//    		args[0] = "-" + option_runTestComputeCellStates + "=test" + File.separator + "quant_006.csv.MarkerStates_Validation.txt";
		

//    		java -jar MOHAtool.jar -runTestHeterogeneity=test/AGA_260_3_AKT_Validation.txt
//    		args = new String[1];
//    		args[0] = "-" + option_runTestHeterogeneity + "=test" + File.separator + "AGA_260_3_AKT_Validation.txt";
    		

//    		java -jar MOHAtool.jar -runTestGeneSetHeterogeneity=test/crcStudy_Validation.txt
    		args = new String[1];
    		args[0] = "-" + option_runTestGeneSetHeterogeneity + "=test" + File.separator + "crcStudy_Validation.txt";


    	}
	    
		// No arguments then default to presenting tool help options
        if (args == null || args.length == 0) {
            args = new String[1];
            args[0] = "-" + option_help;
         }
        
        
    	CommandLineParser commandLine = new CommandLineParser(args);

		verbose = commandLine.getOptionState(option_verbose, false);

		if (commandLine.hasOption(option_h) || commandLine.hasOption(option_help) || commandLine.hasOption(option_questionMark)) {
			printHelp();
		} else if (commandLine.hasOption(option_version)) {
			printVersion();
		} else {
			try {
				initTool(commandLine);
				performToolRun(commandLine);
				closeTool(commandLine);
			} catch (Exception e) {
				System.err.println(e.getMessage());
			}
		}
	}


    
	/**
	 * Class to manage command line parsing and retrieving tool run options and parameter values. 
	 */
    protected class CommandLineParser {
    	
    	protected Map<String, String> cmdOptions = new HashMap<String, String>();
        
    	public CommandLineParser(String[] args) {
    		appendArguments(args);
        }
    	
    	public void appendArguments(String[] args) {
    		if (args != null) {
    			for (int i = 0; i < args.length; i++) {
    				if (args[i].startsWith("-")) {
    					int delimIndx = args[i].indexOf("=");
    					if (delimIndx != -1) {
    						addCommandOption(args[i].substring(1, delimIndx), args[i].substring(delimIndx + 1), true);
    					} else {
    						addCommandOption(args[i].substring(1), "", true);
    					}
    				} else {
    					System.err.println("Command line problem: Do not understand element: " + args[i]);
    				}
    			}
    		}

    	}

    	
        protected void addCommandOption(String optionKey, String optionValue, boolean overwriteExistingOption) {
    		optionKey = optionKey.toLowerCase().trim();
    		optionValue = optionValue.trim();
    		if (cmdOptions.containsKey(optionKey)) {
    			if (overwriteExistingOption) {
    				cmdOptions.remove(optionKey);
    				cmdOptions.put(optionKey, optionValue);
    			}
    		} else {
    			cmdOptions.put(optionKey, optionValue);
    		}
    	}

        public boolean hasOption(String optionKey) {
    		optionKey = optionKey.toLowerCase().trim();
    		return cmdOptions.containsKey(optionKey);
    	}
        
        protected String getOptionValue(String optionKey) {
    		optionKey = optionKey.toLowerCase().trim();
    		return cmdOptions.get(optionKey);
    	}

        public String getOptionValueRemoveQuotes(String optionKey) {
    		String value = getOptionValue(optionKey);
    		if (value != null)
    			return value.replace("\"", "").replace("'", "");
    		else
    			return null;
    	}

        public boolean getOptionState(String optionKey, boolean defaultState) {
    		if (hasOption(optionKey)) {
    			String optionValue = getOptionValue(optionKey);
    			if ((optionValue == null) || (optionValue == "")) {
    				return true;
    			} else if (("F".equalsIgnoreCase(optionValue) || "FALSE".equalsIgnoreCase(optionValue)))
    				return false;
    			else if (("T".equalsIgnoreCase(optionValue) || "TRUE".equalsIgnoreCase(optionValue)))
    				return true;
    			else
    				return defaultState;
    		} else {
    			return defaultState;
    		}
    	}

        public double getOptionDoubleValue(String optionKey, double defaultValue) {
    		if (hasOption(optionKey)) {
    			try {
    				return Double.parseDouble(getOptionValue(optionKey).trim());
    			} catch (Exception e) {
    				System.err.println("Problem parsing option: " + optionKey);
    				return defaultValue;
    			}
    		} else {
    			return defaultValue;
    		}
    	}

        public int getOptionInitValue(String optionKey, int defaultValue) {
    		if (hasOption(optionKey)) {
    			try {
    				return Integer.parseInt(getOptionValue(optionKey).trim());
    			} catch (Exception e) {
    				System.err.println("Problem parsing option: " + optionKey);
    				return defaultValue;
    			}
    		} else {
    			return defaultValue;
    		}
    	}

    }
    
	/**
	 * Parse data line for the given delimiter including special characters. 
	 */
    final protected String[] valueOf(String value, String delimiter) {
		if (delimiter.equals("|"))
			delimiter = "\\|";
		else if (delimiter.equals("."))
			delimiter = "\\.";
		else if (delimiter.equals("("))
			delimiter = "\\(";
		else if (delimiter.equals(")")) 
			delimiter = "\\)";
		return value.split(delimiter);
	}

	/**
	 * Parse data line for the given delimiter and remove double quotations (cvs files) and then any leading or trailing spaces. 
	 */
    final protected String[] valueOfWithTrimAfterRemoveQuotations(String value, String delimiter) {
		String[] values = valueOf(value, delimiter);
		for (int i = 0; i < values.length; i++)
			values[i] = values[i].replace("\"", "").trim();
		return values;
	}
    
	/**
	 * Print out software tool version.
	 */
    protected void printVersion() {
		System.out.println();
		System.out.println("Version: " + this.getClass().getSimpleName() + " " + versionNumber);
		System.out.println();
	}
		
	/**
	 * Print out command line tool usage and options.
	 */
	protected void printHelp() {
		StudyDataTableModel sdtm = new StudyDataTableModel();
		CellDataTableModel cdtm = new CellDataTableModel();

        System.out.println();
        System.out.println("Usage: " + "java -jar MOHAtool.jar " + " [-options]" + " -" + option_computeThresholds  + "=<input Study or Cell Measurement file>" + " -" + option_thresholdFile + "=<output Threshold file>");
        System.out.println("  or");
        System.out.println("       " + "java -jar MOHAtool.jar " + " [-options]" + " -" + option_computeCellStates  + "=<input Study or Cell Measurement file>" + " -" + option_thresholdFile + "=<input Threshold file>" + " -" + option_cmsFile + "=<output Cell Marker State file>" + " -" + option_miFile + "=<output Marker Index file>");
        System.out.println("  or");
        System.out.println("       " + "java -jar MOHAtool.jar " + " [-options]" + " -" + option_computeGeneSetHeterogeneity + "=<input Gene Set file>" + " -" + option_cmsFile + "=<input Cell Marker State file>"  + " -" + option_miFile + "=<input Marker Index file>"  + " -" + option_outputFile + "=<output file>");
        System.out.println("  or");
        System.out.println("       " + "java -jar MOHAtool.jar " + " [-options]" + " -" + option_computeHeterogeneity + "=<input Cell Marker State file>" + " -" + option_outputFile + "=<output file>");
        System.out.println();
        System.out.println("where [-options] command line format is -<option name>=<option value> with a space between each option name-value pair");
        System.out.println();
        System.out.println("where option names include:");
        System.out.println("  -" + option_help + " -" + option_h + " -" + option_questionMark + "       print this help message");
        System.out.println("  -" + option_version + "          print product version and exit");
        System.out.println("  -" + option_verbose + "=false    console output");
        System.out.println("  -" + option_outputDir + "        output directory");
        System.out.println("  -" + option_append + "=true      append output to existing ouput file");
        System.out.println("  -" + option_cdf + "="  + DEFAULT_CRIT_DISTANCE_FACTOR + "         critical distance factor");
        System.out.println("  -" + option_nStateModel + "="  + Integer.toString(DEFAULT_N_STATE_THRESHOLD_MODEL) + "    n-State Marker Threshold Model");
        System.out.println("  -" + option_SampleID + "         sample ID to use if not provided by input Cell file");
        System.out.println("  -" + option_GeneSetName + "      gene (marker) set name to use if not provided by input Gene Set file");
        System.out.println("  -" + option_MaxNumCellStates + " user defined maximum number of molecular cell states");
        System.out.println("  -" + option_MaxNumCells + "=" + cdtm.MAX_NUM_CELLS + " maximum number of cells (memory limited)");
        System.out.println();
        System.out.println("and where Study file default column names are:");
        System.out.println("  -" + option_colName_Sample_ID + "=\"" + sdtm.colName_Sample_ID + "\"" + "          sample ID for cell measure file");
        System.out.println("  -" + option_colName_DATA_FILENAME + "=\"" + sdtm.colName_DATA_FILENAME + "\"" + "  cell measure filename and path");
        System.out.println();
        System.out.println("and where Cell Measurement and Cell Marker State file default column names are:");
        System.out.println("  -" + option_colName_Sample_ID + "=\"" + cdtm.colName_Sample_ID + "\"" + "          optional column");
        System.out.println("  -" + option_colName_Slide_ID + "=\"" + cdtm.colName_Slide_ID + "\"" + "            optional column");
        System.out.println("  -" + option_colName_Position_ID + "=\"" + cdtm.colName_Position_ID + "\"" + "      optional column");
        System.out.println("  -" + option_colName_Cell_ID + "=\"" + cdtm.colName_Cell_ID + "\"" + "              optional column");
        System.out.println("  -" + option_colName_Cell_Center_X + "=\"" + cdtm.colName_Cell_Center_X + "\"" + "  cell spatial X coordinate");
        System.out.println("  -" + option_colName_Cell_Center_Y + "=\"" + cdtm.colName_Cell_Center_Y + "\"" + "  cell spatial Y coordinate");
        System.out.println("  -" + option_colName_Cell_Area + "=\"" + cdtm.colName_Cell_Area + "\"" + "          cell area");
        System.out.println("  -" + option_colName_Cell_Radius + "=\"" + cdtm.colName_Cell_Radius + "\"" + "      cell radius (alternative to cell area)");
        System.out.println("  -" + option_colName_Marker_States + "=\"" + cdtm.colName_Marker_States + "\"" + "  marker states");
        System.out.println("  -" + option_biomarkerColNameTag + "=\"" + cdtm.DEFAULT_BIOMARKER_COLNAME_TAG + "\"" +  "     cell biomarker measurement columns contain this tag");
        
        
	}


	/**
	 * Inititialize tool method called prior to calling performToolRun.
	 * Override this method to perform custom tool setup steps (e.g. open database).
	 */
	protected void initTool(CommandLineParser commandLine) throws Exception {
		// do initialize tool
	}

	/**
	 * Close tool method called after calling performToolRun.
	 * Override this method to perform custom tool shutdown steps (e.g. close database).
	 */
	protected void closeTool(CommandLineParser commandLine) throws Exception {
		// do tool shutdown
	}

	/**
	 * Perform tool run using parameters and options stored within the commandLine object.
	 * 
	 */
	protected void performToolRun(CommandLineParser commandLine) throws Exception {
		long time0 = System.currentTimeMillis();
		
		int nStateThresholdModel = commandLine.getOptionInitValue(option_nStateModel, DEFAULT_N_STATE_THRESHOLD_MODEL);
		if (commandLine.hasOption(option_computeThresholds)) {
			String fileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_computeThresholds);
			CellDataTableModel defaultTable = new CellDataTableModel();
			String biomarkerColNameTag = defaultTable.DEFAULT_BIOMARKER_COLNAME_TAG;
			if (commandLine.hasOption(option_biomarkerColNameTag)) {
				biomarkerColNameTag = commandLine.getOptionValueRemoveQuotes(option_biomarkerColNameTag);
			}
			int maxNumCells = commandLine.getOptionInitValue(option_MaxNumCells, defaultTable.MAX_NUM_CELLS);
			String thresholdfileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_thresholdFile);
			if (thresholdfileNameAndPath == null || thresholdfileNameAndPath.length() == 0) {
				String fn = fileNameAndPath;
				if (fn.endsWith(".txt"))
					fn = fn.replace(".txt", "");
				if (fn.endsWith(".cvs"))
					fn = fn.replace(".cvs", "");
				thresholdfileNameAndPath = fn + ".thresholds.txt";
			}
			performComputeThresholds(fileNameAndPath, nStateThresholdModel, biomarkerColNameTag, maxNumCells, thresholdfileNameAndPath, commandLine);
		} else if (commandLine.hasOption(option_computeCellStates)) {
			String fileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_computeCellStates);
			String thresholdfileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_thresholdFile);
			if (thresholdfileNameAndPath == null || thresholdfileNameAndPath.length() == 0) {
				String fn = fileNameAndPath;
				if (fn.endsWith(".txt"))
					fn = fn.replace(".txt", "");
				if (fn.endsWith(".cvs"))
					fn = fn.replace(".cvs", "");
				thresholdfileNameAndPath = fn + ".thresholds.txt";
			}
			String cmsfileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_cmsFile);
			if (cmsfileNameAndPath == null || cmsfileNameAndPath.trim().length() == 0) {
				cmsfileNameAndPath = fileNameAndPath.replace(".txt", "") + ".MarkerStates.txt";
			}
			String mifileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_miFile);
			if (mifileNameAndPath == null || mifileNameAndPath.trim().length() == 0) {
				mifileNameAndPath = fileNameAndPath.replace(".txt", "") + ".MarkerIndex.txt";
			}
			performComputeCellStates(fileNameAndPath, nStateThresholdModel, thresholdfileNameAndPath, cmsfileNameAndPath, mifileNameAndPath, commandLine);
		} else if (commandLine.hasOption(option_computeGeneSetHeterogeneity)) {
			String gsfileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_computeGeneSetHeterogeneity);
			double criticalDistanceFactor = commandLine.getOptionDoubleValue(option_cdf, DEFAULT_CRIT_DISTANCE_FACTOR);
			String cmsfileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_cmsFile);
			if (cmsfileNameAndPath == null || cmsfileNameAndPath.length() == 0) {
	    		throw new Exception("Provide cell marker state file name and path");
			}
			String mifileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_miFile);
			if (mifileNameAndPath == null || mifileNameAndPath.length() == 0) {
	    		throw new Exception("Provide marker index file name and path");
			}
			String outputFileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_outputFile);
			if (outputFileNameAndPath == null) {
				outputFileNameAndPath = "moha_output";
			}
			String outputDir = commandLine.getOptionValueRemoveQuotes(option_outputDir);
			if (outputDir != null) {
				outputFileNameAndPath = outputDir + File.separator + outputFileNameAndPath;
			}
			boolean appendOutput = commandLine.getOptionState(option_append, true);
			performComputeGeneSetHeterogeneity(gsfileNameAndPath, nStateThresholdModel, criticalDistanceFactor,cmsfileNameAndPath,mifileNameAndPath, outputFileNameAndPath, appendOutput, commandLine);
		} else if (commandLine.hasOption(option_computeHeterogeneity)) {
			String cmsfileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_computeHeterogeneity);
	        String outputFileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_outputFile);
	        if (outputFileNameAndPath == null) {
	        	outputFileNameAndPath = "moha_output";
	        }
	    	String outputDir = commandLine.getOptionValueRemoveQuotes(option_outputDir);
	    	if (outputDir != null) {
	    		outputFileNameAndPath = outputDir + File.separator + outputFileNameAndPath;
	    	}
	        boolean appendOutput = commandLine.getOptionState(option_append, true);
			double criticalDistanceFactor = commandLine.getOptionDoubleValue(option_cdf, DEFAULT_CRIT_DISTANCE_FACTOR);
			performComputeHeterogeneity(cmsfileNameAndPath, nStateThresholdModel,criticalDistanceFactor, outputFileNameAndPath, appendOutput, commandLine);
		} else if (commandLine.hasOption(option_runTestHeterogeneity)) {
			performTestHeterogeneity(commandLine);
		} else if (commandLine.hasOption(option_runTestGeneSetHeterogeneity)) {
			performTestGeneSetHeterogeneity(commandLine);
		} else if (commandLine.hasOption(option_runTestComputeThresholds)) {
			performToolTestComputeThresholds(commandLine);
		} else if (commandLine.hasOption(option_runTestComputeCellStates)) {
			performToolTestComputeCellStates(commandLine);
		} else {
			System.err.println("Did not recognized command");
			printHelp();
		}
		long timeDif = System.currentTimeMillis() - time0;
		System.out.println("Finished tool run (" + (timeDif / 1000) + " seconds)");
	}
	

	/**
	 * Computes thresholds for cell measures.
	 * @param fileNameAndPath           Filename and path string to either a cell measure file or a study file that contains a list of cell measure filenames.
	 * @param nStateThresholdModel      The n-state threshold model to use when creating the threshold.
	 * @param biomarkerColNameTag       biomarkerColNameTag
	 * @param MAX_NUM_CELLS             MAX_NUM_CELLS
	 * @param thresholdfileNameAndPath  thresholdfileNameAndPath
	 * @param commandLine               CommandLineParser object
	 */
	protected void performComputeThresholds(String fileNameAndPath, int nStateThresholdModel, String biomarkerColNameTag, int MAX_NUM_CELLS, String thresholdfileNameAndPath, CommandLineParser commandLine) throws Exception {

		System.out.println("Computing thresholds...");
		
		List<File> cmFiles = getCellMeasureFileList(fileNameAndPath, true, commandLine);
		

		List<String> columnNameList = getColumnNameList(cmFiles);

    	int numCellMeasures = getNumCellsInStudy(cmFiles);
    	

		if (numCellMeasures > MAX_NUM_CELLS) {
			throw new Exception("Number of cell measure for entire study is " + numCellMeasures + " and exceeds current memory limit of " + MAX_NUM_CELLS);
		}


		
		CellMeasurementColumn measurementColumnTag = new CellMeasurementColumn(biomarkerColNameTag);
		
		
		List<CellMeasurementColumn> measurementColNameTags =  measurementColumnTag.generateMeasurementColNameTags(columnNameList, true);

		
		List<CellMeasurementColumn> measurementColumns = new ArrayList<CellMeasurementColumn>();
		for(String columnName : columnNameList) {
			for(CellMeasurementColumn mcTag : measurementColNameTags) {
				if (columnName.toLowerCase().contains(mcTag.measureColNameTag.toLowerCase())) {
					CellMeasurementColumn measurementColumn = new CellMeasurementColumn(mcTag);
					measurementColumn.measurementColName = columnName;
					measurementColumn.biomarkerName = columnName.replace(mcTag.measureColNameTag, "");
					measurementColumns.add(measurementColumn);
				}
			}
		}
		

		List<Threshold> biomarkerThresholds = new ArrayList<Threshold>();
		for(CellMeasurementColumn mc : measurementColumns) {
			
			Threshold threshold = createMarkerThreshold(cmFiles, mc.measurementColName, mc.biomarkerName, mc.location, mc.metric, nStateThresholdModel, numCellMeasures, MAX_NUM_CELLS);
			biomarkerThresholds.add(threshold);
			
			System.out.println("biomarkerColName" + "\t" + threshold.biomarkerColName);
			System.out.println("minValue" + "\t" + threshold.minValue);
			for(int j=0;j<threshold.thresholdValues.length;j++) {
				System.out.println("Threshold_" + (j + 1) + "\t" + threshold.thresholdValues[j]);
			}
			System.out.println("maxValue" + "\t" + threshold.maxValue);
			System.out.println("---------------------------------------");
			//System.out.println("numMeasures" + "\t" + threshold.numCellMeasures);
		}


		saveThresholds(thresholdfileNameAndPath,  nStateThresholdModel, biomarkerThresholds);
			
	}

	
	/**
	 * The cell measurement column class defines the column tags in a cell measurement file.
	 * Each cell measurement column name is composed of a biomarker, location, and metric tag.
	 * This class manages the format which can change from one file to another.
	 */
    protected class CellMeasurementColumn {
    	   
    	public String measurementColName;
    	public String biomarkerName;
    
    	public String measureColNameTag;
    	public String location;
    	public String metric;

    	CellMeasurementColumn(CellMeasurementColumn mcTag) {
        	measurementColName = mcTag.measurementColName;
        	biomarkerName = mcTag.biomarkerName;
        	measureColNameTag = mcTag.measureColNameTag;
        	location = mcTag.location;
        	metric = mcTag.metric;
    	}
    	
		public CellMeasurementColumn(String tag) {
			this.measureColNameTag = tag;
			
			tag = measureColNameTag.toLowerCase();
		    if (tag.contains("median"))
		    	metric = "median";
		    else if (tag.contains("mean"))
		    	metric = "mean";
		    else
		    	metric = measureColNameTag;
		    
		    if (tag.contains("cell"))
		    	location = "cell";
		    else if (tag.contains("cytosol") || tag.contains("cyt") || tag.contains("cytoplasm"))
		    	location = "cytosol";
		    else if (tag.contains("membrane") || tag.contains("memb") || tag.contains("plasma_membrane"))
		    	location = "membrane";
		    else if (tag.contains("nuclear") || tag.contains("nuc"))
		    	location = "nuclear";
		    else if (tag.contains("extracellular") || tag.contains("extracellular_space"))
		    	location = "extracellular";
		    else
		    	location = measureColNameTag;
		}

		public List<String> getLocationStrs() {
			List<String> locations = new ArrayList<String>();
    		locations.add("cell");
    		locations.add("cytosol");
    		locations.add("cytoplasm");
    		locations.add("cyt");
    		locations.add("nuclear");
    		locations.add("nuc");
    		locations.add("membrane");
    		locations.add("memb");
    		locations.add("plasma_membrane");
    		locations.add("extracellular_space");
    		locations.add("extracellular");
    		return locations;
		}
		
		public List<CellMeasurementColumn> generateMeasurementColNameTags(List<String> columnNameList, boolean generateAllLocations) {

			List<CellMeasurementColumn> measurementColNameTags = new ArrayList<CellMeasurementColumn>();
			if (generateAllLocations) {
				
				Set<String> tagSet = new HashSet<String>();
			    String searchMetricTag = "";
			    int indxMetricTagIndx = measureColNameTag.toLowerCase().indexOf(metric);
			    int indxLocationTagIndx = measureColNameTag.toLowerCase().indexOf(location);
			    String ts1;
			    String ts2;
			    String ts3;
			    boolean locationFirst = true;
			    if (indxLocationTagIndx < indxMetricTagIndx) {
			    	locationFirst = true;
				    ts1 = measureColNameTag.substring(0, indxLocationTagIndx);
				    ts2 = measureColNameTag.substring(indxLocationTagIndx + location.length(), indxMetricTagIndx);
				    ts3 = measureColNameTag.substring(indxMetricTagIndx + metric.length(), measureColNameTag.length());
				    searchMetricTag = ts2 + metric + ts3; 

			    } else {
			    	locationFirst = false;
				    ts1 = measureColNameTag.substring(0, indxMetricTagIndx);
				    ts2 = measureColNameTag.substring(indxMetricTagIndx + metric.length(), indxLocationTagIndx);
				    ts3 = measureColNameTag.substring(indxLocationTagIndx + location.length(), measureColNameTag.length());
				    searchMetricTag = ts1 + metric + ts2;
			    }
			    for(String columnName : columnNameList) {
			    	if (columnName.toLowerCase().contains(searchMetricTag.toLowerCase())) {
			    		String searchTag;
			    		for(String locationTag : getLocationStrs()) {
			    		if (locationFirst) {
				    			searchTag = ts1 + locationTag + searchMetricTag; 
				    		} else {
				    			searchTag = searchMetricTag + locationTag + ts3; 
				    		}
				    		if (columnName.toLowerCase().contains(searchTag.toLowerCase())) {
				    			String mcTag = columnName.substring(columnName.length() - searchTag.length(), columnName.length());
				    			tagSet.add(mcTag);
				    			break;
				    		}
			    		}
			    	}
			    }

			    for(String mcTag : tagSet.toArray(new String[0])) {
	    			CellMeasurementColumn tag = new CellMeasurementColumn(mcTag);
	    			measurementColNameTags.add(tag);
			    }
    			
			} else {
				measurementColNameTags.add(this);
			}
			return measurementColNameTags;
		}
		
    }
    

	/**
	 * Threshold class contains the parameters used to convert the continuous measurements into
	 * ordinal values representing distinct n-states for a biomarker at a particular location
	 * and metric type (medium or mean).
	 */
    protected class Threshold {
    	
        public static final String colName_BiomarkerColumnName = "BiomarkerColumnName";
        public static final String colName_BiomarkerName = "BiomarkerName";
        public static final String colName_BiomarkerLocation = "BiomarkerLocation";
    	public static final String colName_BiomarkerMetric = "BiomarkerMetric";
    	public static final String colName_minValue = "minValue";
    	public static final String colName_threshold_ = "threshold_";
    	public static final String colName__of_ = "_of_";
    	public static final String colName_maxValue = "maxValue";
    	
    	public String biomarkerColName;
    	public String biomarkerName;
    	public String biomarkerLocation;
    	public String biomarkerMetric;
    	public double minValue;
    	public double maxValue;
    	public double [] thresholdValues = new double[0];
    	
		public Threshold(int nStateThresholdModel) {
			thresholdValues = new double[nStateThresholdModel-1];
		}

		public int getMarkerStateValue(double markerValue) {
			for(int i = 0;i<thresholdValues.length;i++) {
				if (markerValue < thresholdValues[i]) {
					return i;
				}
			}
			return thresholdValues.length;
		}
    	
    }
     
	/**
	 * An inputed list of Threshold objects is saved to an output file.
	 * @param thresholdfileNameAndPath   The output filename and path to save the list of threshold objects.
	 * @param nStateThresholdModel       The n-state threshold model to load.
	 * @param thresholds                 The list of threshold objects to save to the output file.
	 * @see Threshold
	 */
	protected void saveThresholds(String thresholdfileNameAndPath, int nStateThresholdModel, List<Threshold> thresholds) throws Exception {
		
		boolean append = false;
		PrintWriter pw = new PrintWriter(new FileWriter(thresholdfileNameAndPath, append));
	
	    String dataLine = Threshold.colName_BiomarkerColumnName;
		for(Threshold threshold : thresholds) {
			dataLine += "\t" + threshold.biomarkerColName;
		}
		pw.println(dataLine);
		
	    dataLine = Threshold.colName_BiomarkerName;
		for(Threshold threshold : thresholds) {
			dataLine += "\t" + threshold.biomarkerName;
		}
		pw.println(dataLine);
		
	    dataLine = Threshold.colName_BiomarkerLocation;
		for(Threshold threshold : thresholds) {
			dataLine += "\t" + threshold.biomarkerLocation;
		}
		pw.println(dataLine);
		
	    dataLine = Threshold.colName_BiomarkerMetric;
		for(Threshold threshold : thresholds) {
			dataLine += "\t" + threshold.biomarkerMetric;
		}
		pw.println(dataLine);
		
	    dataLine = Threshold.colName_minValue;
		for(Threshold threshold : thresholds) {
			dataLine += "\t" + threshold.minValue;
		}
		pw.println(dataLine);
		
		for(int i = 0;i<nStateThresholdModel-1;i++) {
		    dataLine = Threshold.colName_threshold_ + Integer.toString(i + 1) + Threshold.colName__of_ + (nStateThresholdModel - 1);
			for(Threshold threshold : thresholds) {
				if (threshold.thresholdValues.length > i)
				dataLine += "\t" + threshold.thresholdValues[i];
				else
					dataLine += "\t" + "NA";	
			}
			pw.println(dataLine);
		}
		
	    dataLine = Threshold.colName_maxValue;
		for(Threshold threshold : thresholds) {
			dataLine += "\t" + threshold.maxValue;
		}
		pw.println(dataLine);
		
		pw.flush();
		pw.close();
		
	}

	/**
	 * Thresholds objects are loaded from an inputted threshold filename and path.
	 * @param thresholdfileNameAndPath  Filename and path to the thresholds to load.
	 * @param nStateThresholdModel      The n-state threshold model to load.
	 * @return Threshold                The List of Threshold objects loaded from the input file.
	 * @see Threshold
	 */
	protected List<Threshold> loadThresholds(String thresholdfileNameAndPath, int nStateThresholdModel) throws Exception {
		
		if (thresholdfileNameAndPath == null || thresholdfileNameAndPath.trim().length() == 0) {
			throw new Exception("Threshold filename and path is not defined.");
		}
		File thresholdFile = new File(thresholdfileNameAndPath);
		if (!thresholdFile.exists()) {
			throw new Exception("Threshold file does not exist: " + thresholdfileNameAndPath);
		}
		
		DefaultTableModel thresholdData = new DefaultTableModel();
    	loadDataTable(thresholdFile, thresholdData, new ArrayList<String>(), new ArrayList<String>(), false);
		if (thresholdData.getRowCount() == 0) {
			throw new Exception("No entries in threshold file: "  + thresholdfileNameAndPath);
		}
		
		if (!Threshold.colName_BiomarkerColumnName.equalsIgnoreCase(thresholdData.getColumnName(0))) {
			throw new Exception("Problem with threshold file. Row 1 column 1 does not equal: " + Threshold.colName_BiomarkerColumnName);
		}
		if (!Threshold.colName_BiomarkerName.equalsIgnoreCase((String)thresholdData.getValueAt(0, 0))) {
			throw new Exception("Problem with threshold file. Row 2 column 1 does not equal: " + Threshold.colName_BiomarkerName);
		}
		if (!Threshold.colName_BiomarkerLocation.equalsIgnoreCase((String)thresholdData.getValueAt(1, 0))) {
			throw new Exception("Problem with threshold file. Row 3 column 1 does not equal: " + Threshold.colName_BiomarkerLocation);
		}
		if (!Threshold.colName_BiomarkerMetric.equalsIgnoreCase((String)thresholdData.getValueAt(2, 0))) {
			throw new Exception("Problem with threshold file. Row 4 column 1 does not equal: " + Threshold.colName_BiomarkerMetric);
		}
		if (!Threshold.colName_minValue.equalsIgnoreCase((String)thresholdData.getValueAt(3, 0))) {
			throw new Exception("Problem with threshold file. Row 5 column 1 does not equal: " + Threshold.colName_minValue);
		}
		for(int i = 0;i<nStateThresholdModel-1;i++) {
			String rowName = Threshold.colName_threshold_ + Integer.toString(i + 1) + Threshold.colName__of_ + (nStateThresholdModel - 1);
			if (!rowName.equalsIgnoreCase((String)thresholdData.getValueAt(4+i, 0))) {
				throw new Exception("Problem with threshold file. Row " +  Integer.toString(6 + i) + " column 1 does not equal: " + rowName);
			}
		}
		if (!Threshold.colName_maxValue.equalsIgnoreCase((String)thresholdData.getValueAt(3 + nStateThresholdModel, 0))) {
			throw new Exception("Problem with threshold file. Row " + Integer.toString(5 + nStateThresholdModel) + " column 1 does not equal: " + Threshold.colName_maxValue);
		}

		List<Threshold> thresholds = new ArrayList<Threshold>();
		for(int ci=1;ci<thresholdData.getColumnCount();ci++) {
			
			Threshold threshold = new Threshold(nStateThresholdModel);
			threshold.biomarkerColName = (String)thresholdData.getColumnName(ci);
			threshold.biomarkerName = (String)thresholdData.getValueAt(0, ci);
			threshold.biomarkerLocation = (String)thresholdData.getValueAt(1, ci);
			threshold.biomarkerMetric = (String)thresholdData.getValueAt(2, ci);
			threshold.minValue = Double.parseDouble((String)thresholdData.getValueAt(3, ci));
			for(int i = 0;i<nStateThresholdModel-1;i++) {
				threshold.thresholdValues[i] = Double.parseDouble((String)thresholdData.getValueAt(4+i, ci));
			}
			threshold.maxValue = Double.parseDouble((String)thresholdData.getValueAt(3 + nStateThresholdModel, ci));

			thresholds.add(threshold);
		}
		
		return thresholds;
		
	}
	
	/**
	 * Returns the total number of cells in inputted list of cell measure files
	 * @param cmFiles   The inputted list of cell measure files.
	 */
    protected int getNumCellsInStudy(List<File> cmFiles) throws Exception {
    	
    	int numCellMeasures = 0;
    	for(File cmFile : cmFiles) {
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(cmFile));
				String dataLine = reader.readLine();
				if (dataLine == null || dataLine.length() < 1) {
					reader.close();
					throw new Exception("First line is blank for cell measure file: " + cmFile.getAbsolutePath());
				}
				while ((dataLine = reader.readLine()) != null) {
					numCellMeasures++;
				}
				reader.close();
			} catch (Exception e) {
				try {
					reader.close();
				} catch (Exception e2) {
				}
				throw e;
			}
		}
    	
    	return numCellMeasures;
    
    } 
    
	/**
	 * Creates a threshold object based on the inputted cell measure files.
	 * @param cmFiles             The inputted list of cell measure files.
	 * @param markerColName          Name of column in the cell measure files to use to create the threshold.
	 * @param markerName             Name of biomarker for the measurement column.
	 * @param markerMetric           The measurement metric (mean, median) for the measurement column.
	 * @param nStateThresholdModel   The n-state threshold model to use when creating the threshold.
	 * @param expectedNumCellMeaures The expected number of cells that will be processed to build the threshold.
	 * @param MAX_NUM_CELLS          The maximum number of cells to process to build the threshold.
	 * @return Threshold             Threshold ceated by the method
	 * @see Threshold
	 */
    protected Threshold createMarkerThreshold(List<File> cmFiles, String markerColName, String markerName, String markerLocation, String markerMetric, int nStateThresholdModel, int expectedNumCellMeaures, int MAX_NUM_CELLS) throws Exception {

    	List<String> requiredColumnNames = new ArrayList<String>();
    	requiredColumnNames.add(markerColName);

    	
		List<Float> data = new ArrayList<Float>(expectedNumCellMeaures);			
				
		int numCellMeasures = 0;
    	for(File file : cmFiles) {
    		
			BufferedReader reader = null;
			try {	
				reader = new BufferedReader(new FileReader(file));
				String dataLine = reader.readLine();
	            //discover column delimiter 
				String delimiter = "\t";
				if (!dataLine.contains("\t"))
					delimiter = ",";
				
	            String[] columnNames = valueOfWithTrimAfterRemoveQuotations(dataLine, delimiter);
				Map<String, Integer> columnMap = new HashMap<String, Integer>();
				for (String colName : requiredColumnNames) {
					boolean found = false;
					for (int i = 0; i < columnNames.length; i++) {
						if (colName.equalsIgnoreCase(columnNames[i])) {
							columnMap.put(colName, i);
							found = true;
							break;
						}
					}
					if (!found) {
						throw new Exception("Required column: " + colName + " was not found in file: " + file.getAbsolutePath());
					}
				}
	

				columnNames = columnMap.keySet().toArray(new String[0]);

				int markerValueColIndx = columnMap.get(markerColName);

				while ((dataLine = reader.readLine()) != null) {
					String[] rowValues = dataLine.split(delimiter);
					try {
						//remove any double quotes
						String value = rowValues[markerValueColIndx].replace("\"", "").trim();
						data.add(Float.parseFloat(value));
						numCellMeasures++;
						if (numCellMeasures >= MAX_NUM_CELLS)
							break;
								
					} catch (Exception e3) {
						// ignore cell measure for this marker if values cannot be parsed to float. For example value = "NA", "na", or "".
					}

				}
	
				reader.close();
				
			} catch (Exception e) {
				try {
					reader.close();
				} catch (Exception e2) {
				}
				throw e;
			}
			
			if (numCellMeasures >= MAX_NUM_CELLS)
				break;
    	}
		

    	Collections.sort(data);
			

		if (data.size() > nStateThresholdModel) {

			Threshold threshold = new Threshold(nStateThresholdModel);
			
			threshold.biomarkerColName = markerColName;
			threshold.biomarkerName = markerName;
			threshold.biomarkerLocation = markerLocation;
			threshold.biomarkerMetric = markerMetric;
			threshold.minValue = data.get(0);
			threshold.maxValue = data.get(data.size() - 1);

			int nzi = 0;
			int m = data.size() / nStateThresholdModel;

			if (Math.abs(data.get(m) - threshold.minValue) < 1e-8f) {
				// value at m is equal to minimum value (e.g. 0 or log2 value)
				// put them all into the lowest threshold bin
				for (nzi = m; nzi < data.size(); nzi++) {
					if ((data.get(nzi) - threshold.minValue) > 1e-8f) {
						break;
					}
				}
				// the first nzv values (get(0) to get(nzi-1) will go into
				// lowest threshold bin
				m = (data.size() - nzi) / nStateThresholdModel;
			}

			
			for (int i = 1; i < nStateThresholdModel; i++) {
				int j = nzi + i * m;
				if (j < data.size()) {
					if (data.size() % nStateThresholdModel == 0) {
						threshold.thresholdValues[i-1] = (data.get(j) + data.get(j - 1)) / 2;
					} else {
						threshold.thresholdValues[i-1] = data.get(j);
					}
				} else {
					threshold.thresholdValues[i-1] = threshold.maxValue;
				}
			}
			
			return threshold;
			
		} else {
			throw new Exception("Not enough cell measures to compute thresholds");
		}
    }
	

	
	/**
	 * Creates a list of column names from a list of files.
	 * @param files   List of files.
	 */
	protected List<String> getColumnNameList(List<File> files) throws Exception {
		
		Set<String> colNameSet = new HashSet<String>();
		List<String> columnNameList = new ArrayList<String>();
		for(File file : files) {
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(file));
				String dataLine = reader.readLine();
				String delimiter = "\t";
				if (!dataLine.contains("\t"))
					delimiter = ",";
	            String[] columnNames = valueOfWithTrimAfterRemoveQuotations(dataLine, delimiter);
				for(String colName : columnNames) {
					if (!colNameSet.contains(colName)) {
						colNameSet.add(colName);
						columnNameList.add(colName);
					}
					
				}
				
				reader.close();
			} catch (Exception e) {
				try {
					reader.close();
				} catch (Exception e2) {
				}
				throw e;
			}
		}

	  return columnNameList;
	}

	
	/**
	 * Class for a study data table model
	 */
	protected class StudyDataTableModel extends DefaultTableModel {
		
		private static final long serialVersionUID = 1L;


		public String colName_Sample_ID = "SAMPLE_ID";
	    public String colName_Slide_ID = "Slide_ID";
	    public String colName_Position_ID = "Position_ID";
		public String colName_PatientID = "PATIENT_ID";
		public String colName_TissueID = "TISSUE_ID";
	    public String colName_DATA_FILENAME = "DATA_FILENAME";
	    
	    public int colIndx_Sample_ID = -1;
	    public int colIndx_Slide_ID = -1;
	    public int colIndx_Position_ID = -1;
		public int colIndx_PatientID= -1;
		public int colIndx_TissueID = -1;
	    public int colIndx_DataFilenameAndPath = -1;

		public StudyDataTableModel() {
			super();
		}
		
    	public List<String> getRequiredColNameList() {
        	List<String> requiredColNameList = new ArrayList<String>();
        	requiredColNameList.add(colName_Sample_ID);
        	requiredColNameList.add(colName_DATA_FILENAME);
    		return requiredColNameList;
    	}
    	
    	public List<String> getAltColNameList() {
        	List<String> altColNameList = new ArrayList<String>();
        	altColNameList.add(colName_Slide_ID);
        	altColNameList.add(colName_Position_ID);
        	altColNameList.add(colName_PatientID);
        	altColNameList.add(colName_TissueID);
    		return altColNameList;
    	}
    	
		public void updateColumnIndexes() {
			colIndx_Sample_ID = findColumn(colName_Sample_ID);
			colIndx_Slide_ID = findColumn(colName_Slide_ID);
			colIndx_Position_ID = findColumn(colName_Position_ID);
			colIndx_PatientID = findColumn(colName_PatientID);
			colIndx_TissueID = findColumn(colName_TissueID);
			colIndx_DataFilenameAndPath = findColumn(colName_DATA_FILENAME);

		}
		
		
		public String getSample_ID(int i) {
			if (colIndx_Sample_ID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_Sample_ID);
		}
		
		public String getSlideID(int i) {
			if (colIndx_Slide_ID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_Slide_ID);
		}
		
		public String getPositionID(int i) {
			if (colIndx_Position_ID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_Position_ID);
		}
		
		public String getPatientID(int i) {
			if (colIndx_PatientID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_PatientID);
		}
		
		public String getTissueID(int i) {
			if (colIndx_TissueID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_TissueID);
		}
		
		public String getDataFilenameAndPath(int i) {
			if (colIndx_DataFilenameAndPath < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_DataFilenameAndPath);
		}
		
	}

	/**
	 * Returns a list of cell measure files from an inputted filename and path.
	 * @param fileNameAndPath
	 * @param assumeCellMeasureFileOnException
	 * @param commandLine
	 */
	protected List<File> getCellMeasureFileList(String fileNameAndPath, boolean assumeCellMeasureFileOnException, CommandLineParser commandLine) throws Exception  {
		
		if (fileNameAndPath == null || fileNameAndPath.length() == 0) {
    		throw new Exception("Provide study file name and path");
		}
		File studyFile = new File(fileNameAndPath);
		if (!studyFile.exists()) {
			throw new Exception("Study or measurement file does not exist: " + fileNameAndPath);
		}
		
		try {
			
			StudyDataTableModel studyData = new StudyDataTableModel();
			if (commandLine.hasOption(option_colName_Sample_ID))
				studyData.colName_Sample_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Sample_ID);
			if (commandLine.hasOption(option_colName_Slide_ID))
				studyData.colName_Slide_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Slide_ID);
			if (commandLine.hasOption(option_colName_Position_ID))
				studyData.colName_Position_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Position_ID);
			if (commandLine.hasOption(option_colName_DATA_FILENAME))
				studyData.colName_DATA_FILENAME = commandLine.getOptionValueRemoveQuotes(option_colName_DATA_FILENAME);
			loadDataTable(studyFile, studyData, studyData.getRequiredColNameList(), studyData.getAltColNameList(), true);
			if (studyData.getRowCount() == 0) {
				throw new Exception("No entries in Study File: "  + fileNameAndPath);
			}
			studyData.updateColumnIndexes();
			
			List<File> cmFiles = getStudyDataFileList(studyData);
			
			return cmFiles;
			
		} catch (Exception e) {
			if (e.toString().contains("Required column") & assumeCellMeasureFileOnException) {
				System.out.println("Assuming input file is a cell measurement file");
				
				List<File> cmFiles = new ArrayList<File>();
				cmFiles.add(studyFile);
				return cmFiles;
			} else {
				// file contains required StudyTableModel column names but there was another problem
				throw e;
			}	
		}
	}
	

	/**
	  * Returns a list of cell measure files from an inputted study data table model
	 * @param studyData            The inputted study data table model.
	 * @see StudyDataTableModel
	 */
	protected List<File> getStudyDataFileList(StudyDataTableModel studyData) throws Exception {
		List<File> fileList = new ArrayList<File>();
    	Set<String> fileSet = new HashSet<String>();
        for(int i=0;i<studyData.getRowCount();i++) {
        	fileSet.add(studyData.getDataFilenameAndPath(i));
        }
        for(String fn : fileSet.toArray(new String[0])) {
        	File file = new File(fn);
			if (!file.exists()) {
				throw new Exception("Cell measures (DataFilenameAndPath) file does not exist: " + fn);
			}
        	fileList.add(file);
        }
    	return fileList;

    }
	

	/**
	 * Updated table model column headings using the provided command line arguments
	 * @param cellData
	 * @param commandLine
	 */
	protected CellDataTableModel updatedTableModelColumnNames(CellDataTableModel cellData, CommandLineParser commandLine) {
		if (commandLine.hasOption(option_colName_Sample_ID))
			cellData.colName_Sample_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Sample_ID);
		if (commandLine.hasOption(option_colName_Slide_ID))
			cellData.colName_Slide_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Slide_ID);
		if (commandLine.hasOption(option_colName_Position_ID))
			cellData.colName_Position_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Position_ID);
		if (commandLine.hasOption(option_colName_Cell_ID))
			cellData.colName_Cell_ID = commandLine.getOptionValueRemoveQuotes(option_colName_Cell_ID);
		if (commandLine.hasOption(option_colName_Cell_Center_X))
			cellData.colName_Cell_Center_X = commandLine.getOptionValueRemoveQuotes(option_colName_Cell_Center_X);
		if (commandLine.hasOption(option_colName_Cell_Center_Y))
			cellData.colName_Cell_Center_Y = commandLine.getOptionValueRemoveQuotes(option_colName_Cell_Center_Y);
		if (commandLine.hasOption(option_colName_Cell_Area))
			cellData.colName_Cell_Area = commandLine.getOptionValueRemoveQuotes(option_colName_Cell_Area);
		if (commandLine.hasOption(option_colName_Cell_Radius))
			cellData.colName_Cell_Radius = commandLine.getOptionValueRemoveQuotes(option_colName_Cell_Radius);
		if (commandLine.hasOption(option_colName_Marker_States))
			cellData.colName_Marker_States = commandLine.getOptionValueRemoveQuotes(option_colName_Marker_States);
		return cellData;
	}
	

	/**
	  * Returns a column index given an inputted column map and a target column name.
	 */
	protected int getColIndx(Map<String, Integer> columnMap, String colName) {
		Integer colIndx = columnMap.get(colName);
		if (colIndx == null)
			colIndx = -1;
		return colIndx;
	}

	/**
	 * Compute cell states from inputted cell measure files, an existing thresholds.
	 * The cell states are outputed and saved to a cell marker state file. A marker index is generted and saved.
	 * @param fileNameAndPath 
	 * @param nStateThresholdModel
	 * @param thresholdfileNameAndPath
	 * @param cmsfileNameAndPath
	 * @param mifileNameAndPath
	 * @param commandLine
	 */
	protected void performComputeCellStates(String fileNameAndPath, int nStateThresholdModel, String thresholdfileNameAndPath, String cmsfileNameAndPath, String mifileNameAndPath, CommandLineParser commandLine) throws Exception {
		
		System.out.println("Computing cell states...");

		List<File> cmFiles = getCellMeasureFileList(fileNameAndPath, true, commandLine);
		

		File thresholdFile = new File(thresholdfileNameAndPath);
		if (!thresholdFile.exists()) {
			throw new Exception("Threshold file does not exist: " + thresholdfileNameAndPath);
		}
		

		List<Threshold> thresholds = loadThresholds(thresholdfileNameAndPath, nStateThresholdModel);
		
		
		File cmsfile = new File(cmsfileNameAndPath);
		if (cmsfile.exists()) {
			//throw new Exception("Please delete existing Marker State File (cmsfile): " + cmsfileNameAndPath);
		}


		File mifile = new File(mifileNameAndPath);
		if (mifile.exists()) {
			//throw new Exception("Please delete existing Marker Index File (mifile): " + mifileNameAndPath);
		}
		

		boolean append = false;
		PrintWriter pw = new PrintWriter(new FileWriter(cmsfile.getAbsolutePath(), append));
		
		
		CellDataTableModel cmsTable = updatedTableModelColumnNames(new CellDataTableModel(), commandLine);
		String headerStr = cmsTable.colName_Slide_ID 
				+ "\t" + cmsTable.colName_Position_ID 
				+ "\t" + cmsTable.colName_Cell_ID 
				+ "\t" + cmsTable.colName_Cell_Center_X 
				+ "\t" + cmsTable.colName_Cell_Center_Y
				+ "\t" + cmsTable.colName_Cell_Area 
				+ "\t" + cmsTable.colName_Marker_States;
		pw.println(headerStr);
		
		
		CellDataTableModel cmTable = updatedTableModelColumnNames(new CellDataTableModel(), commandLine);
    	for(File file : cmFiles) {
    		
			BufferedReader reader = null;
			try {	
				reader = new BufferedReader(new FileReader(file));
				String dataLine = reader.readLine();
	            //discover column delimiter 
				String delimiter = "\t";
				if (!dataLine.contains("\t"))
					delimiter = ",";
				
	            String[] columnNames = valueOfWithTrimAfterRemoveQuotations(dataLine, delimiter);


				Map<String, Integer> columnMap = new HashMap<String, Integer>();
				for (int i = 0; i < columnNames.length; i++) {
					columnMap.put(columnNames[i], i);
				}
				
				List<Integer> biomarkerColIndx = new ArrayList<Integer>();
				for (int i = 0; i < thresholds.size(); i++) {
					biomarkerColIndx.add(getColIndx(columnMap, thresholds.get(i).biomarkerColName));
				}
				
				int colIndx_Slide_ID = getColIndx(columnMap, cmTable.colName_Slide_ID);
				int colIndx_Position_ID = getColIndx(columnMap, cmTable.colName_Position_ID);
				int colIndx_Cell_ID = getColIndx(columnMap, cmTable.colName_Cell_ID);
				int colIndx_Cell_Center_X = getColIndx(columnMap, cmTable.colName_Cell_Center_X);
				int colIndx_Cell_Center_Y = getColIndx(columnMap, cmTable.colName_Cell_Center_Y);
				int colIndx_Cell_Area = getColIndx(columnMap, cmTable.colName_Cell_Area);


				while ((dataLine = reader.readLine()) != null) {
					String[] rowValues = valueOfWithTrimAfterRemoveQuotations(dataLine, delimiter);
					try {

						String slideID = "NA";
						if (colIndx_Slide_ID >= 0) {
							slideID = rowValues[colIndx_Slide_ID];
						}
						
						String spotID = "NA";
						if (colIndx_Position_ID >= 0) {
							spotID = rowValues[colIndx_Position_ID];
						}
						
						String cellID = "NA";
						if (colIndx_Cell_ID >= 0) {
							cellID = rowValues[colIndx_Cell_ID];
						}
						
						String cellCenterX = "NA";
						if (colIndx_Cell_Center_X >= 0) {
							cellCenterX = rowValues[colIndx_Cell_Center_X];
						}
						
						String cellCenterY = "NA";
						if (colIndx_Cell_Center_Y >= 0) {
							cellCenterY = rowValues[colIndx_Cell_Center_Y];
						}
						

						String cellArea;
						if (colIndx_Cell_Area >= 0) {
							cellArea = rowValues[colIndx_Cell_Area];
						} else {
							cellArea = Integer.toString(DEFAULT_CELL_AREA);
						}
						
						
						StringBuffer markerStates = new StringBuffer();
						for (int i = 0; i < thresholds.size(); i++) {
							if (biomarkerColIndx.get(i) >= 0) {
								try {
									// remove any double quotes
									String value = rowValues[biomarkerColIndx.get(i)].replace("\"", "");
									Threshold biomarkerThreshold = thresholds.get(i);
									markerStates.append(biomarkerThreshold.getMarkerStateValue(Double.parseDouble(value)));

								} catch (Exception e3) {
									// ignore cell measure for this marker if
									// values cannot be parsed to float. For
									// example value = "NA", "na", or "".
									markerStates.append("X");
								}
							} else {
								markerStates.append("X");
							}
						}
						

						String dataLineStr = slideID 
								+ "\t" + spotID 
								+ "\t" + cellID 
								+ "\t" + cellCenterX
								+ "\t" + cellCenterY 
								+ "\t" + cellArea 
								+ "\t" + markerStates.toString();
						
						pw.println(dataLineStr);

								
					} catch (Exception e3) {
						// ignore cell measure for this marker if values cannot be parsed to float. For example value = "NA", "na", or "".
					}

				}
	
				reader.close();
				
			} catch (Exception e) {
				try {
					reader.close();
				} catch (Exception e2) {
				}
				throw e;
			}
			
    	}

		pw.flush();
		pw.close();
		
		
		
		pw = new PrintWriter(new FileWriter(mifile.getAbsolutePath(), append));
		
		MarkerIndexTableModel miTable = new MarkerIndexTableModel();
		
	    headerStr = miTable.colName_MarkerIndx 
	    		+ "\t" + miTable.colName_MarkerName 
	    		+ "\t" + miTable.colName_LocationName 
	    		+ "\t" + miTable.colName_MarkerTargets;
	    
		pw.println(headerStr);
		
		for(int i=0;i<thresholds.size();i++) {
			
			// currently setting MarkerTargets column to be the name of the biomarker.  Should be the gene symbols representing the marker.
			String markerTargets = thresholds.get(i).biomarkerName + " "  + thresholds.get(i).biomarkerLocation;

			String dataLineStr = Integer.toString(i) 
					+ "\t" + thresholds.get(i).biomarkerName 
					+ "\t" + thresholds.get(i).biomarkerLocation 
					+ "\t" + markerTargets;
			
			pw.println(dataLineStr);
		}
		
		pw.flush();
		
		pw.close();

	}
	

	/**
	 * Class for a marker index data table model
	 */
	protected class MarkerIndexTableModel extends DefaultTableModel {
		
		private static final long serialVersionUID = 1L;
		
		public String colName_MarkerIndx = "MarkerIndx";
	    public String colName_MarkerName = "MarkerName";
	    public String colName_LocationName = "LocationName";
	    public String colName_MarkerTargets = "MarkerTargets";
	    
	    public int colIndx_MarkerIndx = -1;
	    public int colIndx_MarkerName = -1;
	    public int colIndx_LocationName = -1;
		public int colIndx_MarkerTargets= -1;

		
		public MarkerIndexTableModel() {
			super();
		}
		
    	public List<String> getRequiredColNameList() {
        	List<String> requiredColNameList = new ArrayList<String>();
        	requiredColNameList.add(colName_MarkerName);
        	requiredColNameList.add(colName_MarkerTargets);
    		return requiredColNameList;
    	}
    	
    	public List<String> getAltColNameList() {
        	List<String> altColNameList = new ArrayList<String>();
        	altColNameList.add(colName_MarkerIndx); // if not present assume they are in order found in column msFileColName_Marker_States
        	altColNameList.add(colName_LocationName);
    		return altColNameList;
    	}
    	

		public void updateColumnIndexes() {
			colIndx_MarkerIndx = findColumn(colName_MarkerIndx);
			colIndx_MarkerName = findColumn(colName_MarkerName);
			colIndx_LocationName = findColumn(colName_LocationName);
			colIndx_MarkerTargets = findColumn(colName_MarkerTargets);
		}
		
		public String getMarkerTargets(int i) {
			return (String) getValueAt(i, colIndx_MarkerTargets);
		}
		
	}
	
	

	/**
	 * Class for gene set data table model
	 */
	protected class GeneSetTableModel extends DefaultTableModel {
		
		private static final long serialVersionUID = 1L;
		
		public String geneset_delimiter = ":";
		
		public String colName_markerSetName = "markerSetName";
	    public String colName_markers = "markers";

	    public int colIndx_markerSetName = -1;
	    public int colIndx_markers = -1;

		
		public GeneSetTableModel() {
			super();
		}
		
    	public List<String> getRequiredColNameList() {
        	List<String> requiredColNameList = new ArrayList<String>();
        	requiredColNameList.add(colName_markerSetName);
        	requiredColNameList.add(colName_markers);
    		return requiredColNameList;
    	}
    	
    	public List<String> getAltColNameList() {
        	List<String> altColNameList = new ArrayList<String>();
    		return altColNameList;
    	}
    	

		public void updateColumnIndexes() {
			colIndx_markerSetName = findColumn(colName_markerSetName);
			colIndx_markers = findColumn(colName_markers);
		}
		
		public String getMarkerSetName(int i) {
			return (String) getValueAt(i, colIndx_markerSetName);
		}
		
		public String getMarkers(int i) {
			return (String) getValueAt(i, colIndx_markers);
		}
	}
	

	/**
	 * Computes the heterogeneity metrics for cells for inputted gene sets (i.e. combinations of markers). The cell marker state data is inputted along with a marker index file.
	 * @param gsfileNameAndPath
	 * @param nStateThresholdModel
	 * @param criticalDistanceFactor
	 * @param cmsfileNameAndPath
	 * @param mifileNameAndPath
	 * @param outputFileNameAndPath
	 * @param appendOutput
	 * @param commandLine
	 */
	protected void performComputeGeneSetHeterogeneity(String gsfileNameAndPath, int nStateThresholdModel, double criticalDistanceFactor, String cmsfileNameAndPath, String mifileNameAndPath, String outputFileNameAndPath, boolean appendOutput, CommandLineParser commandLine) throws Exception {
		
		System.out.println("Computing heterogeneity metrics of gene sets...");
		
		if (gsfileNameAndPath == null || gsfileNameAndPath.length() == 0) {
    		throw new Exception("Provide gene set file name and path");
		}
		File gsfile = new File(gsfileNameAndPath);
		if (!gsfile.exists()) {
			throw new Exception("Gene Set File (gsfile) does not exist: " + gsfileNameAndPath);
		}

		File cmsfile = new File(cmsfileNameAndPath);
		if (!cmsfile.exists()) {
			throw new Exception("Cell Marker State File (cmsfile) does not exist: " + cmsfileNameAndPath);
		}

		File mifile = new File(mifileNameAndPath);
		if (!mifile.exists()) {
			throw new Exception("Marker Index File (mifile) does not exist: " + mifileNameAndPath);
		}

    	CellDataTableModel cmsData = updatedTableModelColumnNames(new CellDataTableModel(), commandLine);
    	loadDataTable(cmsfile, cmsData, cmsData.getRequiredColNameList(), cmsData.getAltColNameList(), true);
		if (cmsData.getRowCount() == 0) {
			throw new Exception("No cells (rows) in Cell Marker State File (cmsfile): "  + cmsfileNameAndPath);
		}
		cmsData.parseAndValidateData(cmsfileNameAndPath);
		cmsData.updateColumnIndexes();
		
    	MarkerIndexTableModel miData = new MarkerIndexTableModel();
    	loadDataTable(mifile, miData, miData.getRequiredColNameList(), miData.getAltColNameList(), true);
		if (miData.getRowCount() == 0) {
			throw new Exception("No markers (rows) in Marker Index File (mifile): "  + mifileNameAndPath);
		}
		miData.updateColumnIndexes();

    	GeneSetTableModel gsData = new GeneSetTableModel();
    	loadDataTable(gsfile, gsData, gsData.getRequiredColNameList(), gsData.getAltColNameList(), true);
		if (gsData.getRowCount() == 0) {
			throw new Exception("No Gene Sets (rows) in Gene Set File (gsfile): "  + gsfileNameAndPath);
		}
		gsData.updateColumnIndexes();
		
		
		Set<String> sampleIDset = new HashSet<String>();
    	for(int msi = 0;msi<cmsData.getRowCount();msi++) {
    		sampleIDset.add(cmsData.getSampleID(msi));
    	}
    	List<String> sampleIDs = Arrays.asList(sampleIDset.toArray(new String[0]));
    	Collections.sort(sampleIDs);
		
		
		for (int gsi = 0; gsi < gsData.getRowCount(); gsi++) {
			String geneSetName = gsData.getMarkerSetName(gsi);
			String geneSetMarkerStr = gsData.getMarkers(gsi);

			String [] geneSetMarkerList =  valueOfWithTrimAfterRemoveQuotations(geneSetMarkerStr, gsData.geneset_delimiter);
			for(int i=0;i<geneSetMarkerList.length;i++) {
				if (!geneSetMarkerList[i].contains(" ")) {
					geneSetMarkerList[i] += " cell";
				}
			}

			List<Integer> markerIndexList = new ArrayList<Integer>();
			for(int mii=0;mii<miData.getRowCount();mii++) {
				String markerTargets =  miData.getMarkerTargets(mii);
				String [] msl =  valueOfWithTrimAfterRemoveQuotations(markerTargets, gsData.geneset_delimiter);
				for(int msli=0;msli<msl.length;msli++) {
					for(int si=0;si<geneSetMarkerList.length;si++) {
						 if ((geneSetMarkerList[si].equalsIgnoreCase(msl[msli]))) {
							 markerIndexList.add(mii);
						 }
					}
				}
			}
		
			if (markerIndexList.size() == 0) {
				System.out.println("geneSetName" + "=" + geneSetName);
				System.out.println("numMarkers" + "=" + markerIndexList.size());
				continue;
			}

			for (String sampleID : sampleIDs) {
				
				Vector<String> csCol_Slide_ID = new Vector<String>();
				Vector<String> csCol_Position_ID = new Vector<String>();
				Vector<Double> csCol_Cell_Center_X = new Vector<Double>();
				Vector<Double> csCol_Cell_Center_Y = new Vector<Double>();
				Vector<Double> csCol_Cell_Radius = new Vector<Double>();
				Vector<String> csCol_Marker_States = new Vector<String>();
			    int msi = -1;
		    	try {
			    	for(msi = 0;msi<cmsData.getRowCount();msi++) {
			    		if (sampleID.equals(cmsData.getSampleID(msi))) {
			    			
							String cellMarkerStates = cmsData.getCellMarkerStates(msi);
							char[] markerStatesCharArray = cellMarkerStates.toCharArray();
							StringBuilder sb = new StringBuilder();
							for (int i : markerIndexList) {
								sb.append(markerStatesCharArray[i]);
							}
			    			String cellState = sb.toString();
			    			// do not use cell if cell state contains marker level that is unknown
			    			if (!cellState.contains("x") && !cellState.contains("X")) {
					    		csCol_Slide_ID.addElement(cmsData.getSlideID(msi));
								csCol_Position_ID.addElement(cmsData.getPositionID(msi));
								csCol_Cell_Center_X.addElement(cmsData.getCellCenterX(msi));
								csCol_Cell_Center_Y.addElement(cmsData.getCellCenterY(msi));
								csCol_Cell_Radius.addElement(cmsData.getCell_Radius(msi));
								csCol_Marker_States.addElement(cellState);
			    			}
			    		}
			    	}
		    	} catch (Exception e) {
		    		throw new Exception("Problem with processing cell state data for row " + (msi + 1) + " in file " + cmsfileNameAndPath);
		    	}
		    	CellDataTableModel cellData = updatedTableModelColumnNames(new CellDataTableModel(), commandLine);
				cellData.addColumn(cellData.colName_Slide_ID, csCol_Slide_ID);
				cellData.addColumn(cellData.colName_Position_ID, csCol_Position_ID);
				cellData.addColumn(cellData.colName_Cell_Center_X, csCol_Cell_Center_X);
				cellData.addColumn(cellData.colName_Cell_Center_Y, csCol_Cell_Center_Y);
				cellData.addColumn(cellData.colName_Cell_Radius, csCol_Cell_Radius);
				cellData.addColumn(cellData.colName_Marker_States, csCol_Marker_States);

				int numCells = cellData.getRowCount();

				String mohaContext = geneSetName;

				int maxNumMolecularStates = 1;
				if (commandLine.hasOption(option_MaxNumCellStates)) {
					 maxNumMolecularStates  = commandLine.getOptionInitValue(option_MaxNumCellStates, maxNumMolecularStates);
				} else {
					int numMarkers = markerIndexList.size();
					maxNumMolecularStates = (int)Math.pow(nStateThresholdModel, numMarkers);
				}
				
				System.out.println("sampleID" + "=" + sampleID);
				System.out.println("numCells" + "=" + numCells);
				System.out.println("geneSetName" + "=" + geneSetName);
				System.out.println("numMarkers" + "=" + markerIndexList.size());
				if (verbose) {
					System.out.println("maxNumCellStates" + "=" + maxNumMolecularStates);
					System.out.println(option_cdf + "=" + criticalDistanceFactor);
					System.out.println(option_SampleID + "=" + sampleID);
					System.out.println(option_GeneSetName + "=" + mohaContext);
					System.out.println(option_nStateModel + "=" + nStateThresholdModel);
					System.out.println(option_MaxNumCellStates + "=" + maxNumMolecularStates);
				    System.out.println(option_outputFile + "=" + outputFileNameAndPath);
				    System.out.println(option_append + "=" + appendOutput);
				    
				}
				System.out.println("-------------------------------------------");
				
				List<DiversityRec> diversityRecs = computeMOHAmetrics(cellData, maxNumMolecularStates, criticalDistanceFactor);
				
				printDiversityMetrics(sampleID, mohaContext, diversityRecs, numCells,  maxNumMolecularStates, df4, outputFileNameAndPath, appendOutput);
				
				
			}
		}
		
	}




	/**
	 * Computes the heterogeneity metrics for an inputted cell state file
	 * @param cmsfileNameAndPath
	 * @param nStateThresholdModel
	 * @param criticalDistanceFactor
	 * @param outputFileNameAndPath
	 * @param appendOutput
	 * @param commandLine
	 */
	protected void performComputeHeterogeneity(String cmsfileNameAndPath, int nStateThresholdModel, double criticalDistanceFactor, String outputFileNameAndPath, boolean appendOutput, CommandLineParser commandLine) throws Exception {

		System.out.println("Computing heterogeneity metrics...");
		
		
		if (cmsfileNameAndPath == null || cmsfileNameAndPath.length() == 0) {
    		throw new Exception("Provide cell state file name and path");
		}
		// make sure data table files exist
    	File cmsfile = new File(cmsfileNameAndPath);
    	if (!cmsfile.exists()) {
    		throw new Exception("Cell Marker State File (cmsfile) does not exist: " + cmsfileNameAndPath);
    	}	
		
    	CellDataTableModel cmsData = updatedTableModelColumnNames(new CellDataTableModel(), commandLine);
    	loadDataTable(cmsfile, cmsData, cmsData.getRequiredColNameList(), cmsData.getAltColNameList(), true);
		if (cmsData.getRowCount() == 0) {
			throw new Exception("No cells (rows) in Cell Marker State File (cmsfile): "  + cmsfileNameAndPath);
		}
		cmsData.parseAndValidateData(cmsfileNameAndPath);
		cmsData.updateColumnIndexes();
		
    	
    	int numCells = cmsData.getRowCount();
    	
		String sampleID = commandLine.getOptionValueRemoveQuotes(option_SampleID);
		if (sampleID == null) {
			sampleID = cmsfile.getName();
		}

		String mohaContext = commandLine.getOptionValueRemoveQuotes(option_GeneSetName);
		if (mohaContext == null) {
			mohaContext = cmsfile.getName();
		}

        int maxNumMolecularStates = 1;
        if (commandLine.hasOption(option_MaxNumCellStates)) {
        	 maxNumMolecularStates  = commandLine.getOptionInitValue(option_MaxNumCellStates, maxNumMolecularStates);
        } else {
        	String exampleState = (String)cmsData.getValueAt(0,  cmsData.findColumn(cmsData.colName_Marker_States));
        	int numMarkers = exampleState.length();
        	maxNumMolecularStates = (int)Math.pow(nStateThresholdModel, numMarkers);
        }

    	
    	System.out.println(option_computeHeterogeneity + "=" + cmsfileNameAndPath);
    	System.out.println("numCells" + "=" + numCells);
        if (verbose) {
    		System.out.println("maxNumCellStates" + "=" + maxNumMolecularStates);
    		System.out.println(option_cdf + "=" + criticalDistanceFactor);
    		System.out.println(option_SampleID + "=" + sampleID);
    		System.out.println(option_GeneSetName + "=" + mohaContext);
    		System.out.println(option_nStateModel + "=" + nStateThresholdModel);
    		System.out.println(option_MaxNumCellStates + "=" + maxNumMolecularStates);
            System.out.println(option_outputFile + "=" + outputFileNameAndPath);
            System.out.println(option_append + "=" + appendOutput);
            
    	}
    	
        List<DiversityRec> diversityRecs = computeMOHAmetrics(cmsData, maxNumMolecularStates, criticalDistanceFactor);

        printDiversityMetrics(sampleID, mohaContext, diversityRecs, numCells,  maxNumMolecularStates, df4, outputFileNameAndPath, appendOutput);

	}
	
	/**
	 * Data is loaded by this method into a DefaultTableModel object.
	 * @param file
	 * @param dataTable
	 * @param requiredColumnNames
	 * @param altColumnNames
	 * @param onlyLoadListedColumns
	 * @return DefaultTableModel
	 * @see DefaultTableModel
	 */
    protected DefaultTableModel loadDataTable(File file, DefaultTableModel dataTable, List<String> requiredColumnNames, List<String> altColumnNames, boolean onlyLoadListedColumns) throws Exception {

	
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            String dataLine = reader.readLine();
            
            //discover column delimiter 
			String delimiter = "\t";
			if (!dataLine.contains("\t"))
				delimiter = ",";
			
            String[] columnNames = valueOfWithTrimAfterRemoveQuotations(dataLine, delimiter);

        	int numDataColumns = columnNames.length;
        	
            Map<String,Integer> columnMap = new HashMap<String,Integer>();
            for(String colName : requiredColumnNames) {
            	boolean found = false;
            	for(int i=0;i<columnNames.length;i++) {
            		if (colName.equalsIgnoreCase(columnNames[i])) {
            			columnMap.put(colName, i);
            			found = true;
            			break;
            		}
            	}
            	if (!found) {
            		reader.close();
            		throw new Exception("Required column: " + colName + " was not found in file: " + file.getAbsolutePath());
            	}
            }
            for(String colName : altColumnNames) {
            	for(int i=0;i<columnNames.length;i++) {
            		if (colName.equalsIgnoreCase(columnNames[i])) {
            			columnMap.put(colName, i);
            			break;
            		}
            	}
            }
		
            
            if (onlyLoadListedColumns) {
            	columnNames = columnMap.keySet().toArray(new String[0]);
            	for(String columnName : columnNames) {
            		dataTable.addColumn(columnName);
            	}
            	int [] requiredRowIndexes = new int[columnNames.length];
                for(int i=0;i<columnNames.length;i++) {
                	requiredRowIndexes[i] = columnMap.get(columnNames[i]);
                }
            	String[] requiredRowValues = new String[columnNames.length];
                while ((dataLine = reader.readLine()) != null) {
                    String[] rowValues = dataLine.split(delimiter);
                    if (rowValues.length >= numDataColumns) {
	                    for(int i=0;i<columnNames.length;i++) {
	                    	//remove quoations and trim ends
	                    	requiredRowValues[i] = rowValues[requiredRowIndexes[i]].replace("\"", "").trim();
	                    }
	                    dataTable.addRow(requiredRowValues);
                    }
                }
            } else {
            	for(String columnName : columnNames) {
            		dataTable.addColumn(columnName);
            	}
                while ((dataLine = reader.readLine()) != null) {
                    String[] rowValues = valueOfWithTrimAfterRemoveQuotations(dataLine, delimiter);
                    if (rowValues.length >= numDataColumns) {
	                    dataTable.addRow(rowValues);
                    }
                }	
            }

            reader.close();
            
            return dataTable;
            
        } catch (Exception e) {
            try {
                reader.close();
            } catch (Exception e2) {  }
            throw e;
        }
    }
	
    
	
	private void printTestCommandLine(String[] args) {
  		StringBuffer sb = new StringBuffer();
		sb.append("java -jar MOHAtool.jar ");
		for(String arg : args) {
			sb.append(" ");
			sb.append(arg);
		}
		
		System.out.println("--------------------------------------------");
		System.out.println("Run command line:");
		System.out.println("");
		System.out.println(sb.toString());
		System.out.println("");
		System.out.println("--------------------------------------------");
	}
	
	
	/**
	 * Perform an internal test on computing Heterogeneity metrics.
	 * @param commandLine
	 */
    protected void performTestHeterogeneity(CommandLineParser commandLine) {
		
		System.out.println("--------------------------------------------");
		System.out.println("Running in test mode");
		System.out.println("Compute heterogeneity metrics from cell marker state file");
		System.out.println("--------------------------------------------");
		
    	String validationFileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_runTestHeterogeneity);
    	String csInputFileNameAndPath = validationFileNameAndPath.replace("_Validation.txt", ".MarkerStates.txt");
    	String outputFileNameAndPath = validationFileNameAndPath.replace("_Validation.txt", ".moha.txt");

    	File outputFile = new File(outputFileNameAndPath);
    	if (outputFile.exists()) {
    		outputFile.delete();
    	}
    	
    	String [] argsToAppend = new String[6];
		argsToAppend[0] = "-" + option_computeHeterogeneity +  "=" +  csInputFileNameAndPath;
		argsToAppend[1] = "-" + option_outputFile + "=" + outputFileNameAndPath;
		argsToAppend[2] = "-" + option_append + "=" + "true";
		argsToAppend[3] = "-" + option_SampleID + "=" + "AGA_260_3_222";
		argsToAppend[4] = "-" + option_GeneSetName + "=" + "AKT_Dec2015";
		argsToAppend[5] = "-" + option_cdf + "=" +"1.31";
		
		commandLine.appendArguments(argsToAppend);
		printTestCommandLine(argsToAppend);
		


        try {
            initTool(commandLine);
            performToolRun(commandLine);
            closeTool(commandLine);
            validateToolRun(outputFileNameAndPath, validationFileNameAndPath);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
        
	}
	
	/**
	 * Perform an internal test on computing Heterogeneity metrics on cells for combinations of markers from gene sets.
	 * @param commandLine
	 */
    protected void performTestGeneSetHeterogeneity(CommandLineParser commandLine) {
		
		System.out.println("--------------------------------------------");
		System.out.println("Running in test mode");
		System.out.println("Compute heterogeneity metrics for gene sets from cell marker state file");
		System.out.println("--------------------------------------------");
		
    	String validationFileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_runTestGeneSetHeterogeneity);
    	String gsInputFileNameAndPath = validationFileNameAndPath.replace("_Validation.txt", ".GeneSets.txt");
    	String msInputFileNameAndPath = validationFileNameAndPath.replace("_Validation.txt", ".MarkerStates.txt");
    	String miInputFileNameAndPath = validationFileNameAndPath.replace("_Validation.txt", ".MarkerIndex.txt");
    	String outputFileNameAndPath = validationFileNameAndPath.replace("_Validation.txt", ".moha.txt");
    	
    	File outputFile = new File(outputFileNameAndPath);
    	if (outputFile.exists()) {
    		outputFile.delete();
    	}
    	
    	String [] argsToAppend = new String[6];
		argsToAppend[0] = "-" + option_computeGeneSetHeterogeneity + "=" +  gsInputFileNameAndPath;
		argsToAppend[1] = "-" + option_cmsFile + "=" +  msInputFileNameAndPath;
		argsToAppend[2] = "-" + option_miFile + "=" +  miInputFileNameAndPath;
		argsToAppend[3] = "-" + option_outputFile +"=" + outputFileNameAndPath;
		argsToAppend[4] = "-" + option_append + "=" + "true";
		argsToAppend[5] = "-" + option_cdf + "=" +"1.31";

		commandLine.appendArguments(argsToAppend);
		printTestCommandLine(argsToAppend);
		
        try {
            initTool(commandLine);
            performToolRun(commandLine);
            closeTool(commandLine);
            validateToolRun(outputFileNameAndPath, validationFileNameAndPath);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

	}

	
	/**
     * Perform an internal test on computing thresholds from a cell measurement file.
	 * @param commandLine
	 */
    protected void performToolTestComputeThresholds(CommandLineParser commandLine) {
		
		System.out.println("--------------------------------------------");
		System.out.println("Running in test mode");
		System.out.println("Compute thresholds from cell measure file");
		System.out.println("--------------------------------------------");
		
    	String validationFileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_runTestComputeThresholds);
    	String cmInputFileNameAndPath = validationFileNameAndPath.replace(".thresholds_Validation.txt", "");
    	String outputFileNameAndPath = validationFileNameAndPath.replace(".thresholds_Validation.txt", ".thresholds.txt");

    	File outputFile = new File(outputFileNameAndPath);
    	if (outputFile.exists()) {
    		outputFile.delete();
    	}
    	
    	String [] argsToAppend = new String[2];
		argsToAppend[0] = "-" + option_computeThresholds +  "=" +  cmInputFileNameAndPath;
		argsToAppend[1] = "-" + option_biomarkerColNameTag + "=" + "_Cell_Median";

		commandLine.appendArguments(argsToAppend);
		printTestCommandLine(argsToAppend);

        try {
            initTool(commandLine);
            performToolRun(commandLine);
            closeTool(commandLine);
            validateToolRun(outputFileNameAndPath, validationFileNameAndPath);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

	}
	
	/**
     * Perform an internal test on computing the cell state from an inputted cell measurement file and threshold file.
	 * @param commandLine
	 */
    protected void performToolTestComputeCellStates(CommandLineParser commandLine) {
		System.out.println("--------------------------------------------");
		System.out.println("Running in test mode");
		System.out.println("Compute cell states from cell measure and theshold file");
		System.out.println("--------------------------------------------");
		
    	String validationFileNameAndPath = commandLine.getOptionValueRemoveQuotes(option_runTestComputeCellStates);
    	String cmInputFileNameAndPath = validationFileNameAndPath.replace(".MarkerStates_Validation.txt", "");
       	String thresholdInputFileNameAndPath = validationFileNameAndPath.replace(".MarkerStates_Validation.txt", ".thresholds.txt");
    	String cmsOutputFileNameAndPath = validationFileNameAndPath.replace(".MarkerStates_Validation.txt", ".MarkerStates.txt");
    	String miOutputFileNameAndPath = validationFileNameAndPath.replace(".MarkerStates_Validation.txt", ".MarkerIndex.txt");

    	File cmsOutputFile = new File(cmsOutputFileNameAndPath);
    	if (cmsOutputFile.exists()) {
    		cmsOutputFile.delete();
    	}
    	File  miOutputFile = new File(miOutputFileNameAndPath);
    	if (miOutputFile.exists()) {
    		miOutputFile.delete();
    	}
    	
    	String [] argsToAppend = new String[2];
		argsToAppend[0] = "-" + option_computeCellStates +  "=" +  cmInputFileNameAndPath;
		argsToAppend[1] = "-" + option_thresholdFile + "=" + thresholdInputFileNameAndPath;

		commandLine.appendArguments(argsToAppend);
		printTestCommandLine(argsToAppend);


        try {
            initTool(commandLine);
            performToolRun(commandLine);
            closeTool(commandLine);
            validateToolRun(cmsOutputFileNameAndPath, validationFileNameAndPath);
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
	}
		
	/**
	 * Method compares the contents of two data table files.
	 * @param testFileNameAndPath
	 * @param validationFileNameAndPath
	 */
    protected void validateToolRun(String testFileNameAndPath, String validationFileNameAndPath) throws Exception {
    	
		File testfile = new File(testFileNameAndPath);
    	if (!testfile.exists()) {
    		throw new Exception("Test file does not exist: " + testFileNameAndPath);
    	}	
    	
    	File validationfile = new File(validationFileNameAndPath);
    	if (!validationfile.exists()) {
    		throw new Exception("Validation file does not exist: " + validationFileNameAndPath);
    	}	
    	
    	DefaultTableModel testData = new DefaultTableModel();
    	loadDataTable(testfile, testData, new ArrayList<String>(), new ArrayList<String>(), false);
    	
    	DefaultTableModel validationData = new DefaultTableModel();
    	loadDataTable(validationfile, validationData, new ArrayList<String>(), new ArrayList<String>(), false);
    	
    	
    	boolean problemFound = false;
       	System.out.println();
       	if (testData.getColumnCount() == validationData.getColumnCount()) {
	    	if (testData.getRowCount() == validationData.getRowCount()) {
        		for(int j=0;j<testData.getColumnCount();j++) {
        			if (!testData.getColumnName(j).equals(validationData.getColumnName(j))) {
        				System.out.println("TEST PROBLEM: Column names are not equivalent for column " + (j+1));
        				problemFound = true;
        			}
        			for(int i=0;i<testData.getRowCount();i++) {

	        			if (!testData.getValueAt(i, j).equals(validationData.getValueAt(i, j))) {
	        				System.out.println("TEST PROBLEM: Data values are not equivalent for data row " + (i + 1) + " column " + (j+1));
	        				problemFound = true;
	        			}
	        		}
	    		}
	       	} else {
	       		System.out.println("TEST PROBLEM: Number of data rows are not the same.");
	       		problemFound = true;
	       	}
       	} else {
       		System.out.println("TEST PROBLEM: Number of data columns are not the same.");
       		problemFound = true;
       	}

       	if (problemFound) {
       		System.out.println(this.getClass().getSimpleName() + " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
       		System.out.println(this.getClass().getSimpleName() + " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
       		System.out.println(this.getClass().getSimpleName() + "  TEST COMPLETED - PROBLEMS FOUND   ");
       		System.out.println(this.getClass().getSimpleName() + " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
       		System.out.println(this.getClass().getSimpleName() + " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
       	} else {
       		System.out.println(this.getClass().getSimpleName() + " +++++++++++++++++++++++++++++++++++");
       		System.out.println(this.getClass().getSimpleName() + " +++++++++++++++++++++++++++++++++++");
       		System.out.println(this.getClass().getSimpleName() + "     TEST COMPLETED - SUCCESS!      ");
       		System.out.println(this.getClass().getSimpleName() + " +++++++++++++++++++++++++++++++++++");
       		System.out.println(this.getClass().getSimpleName() + " +++++++++++++++++++++++++++++++++++");
       	}
	}
		
    // ===============================================================================================================================
	// Multi-Omics Heterogeneity Analysis (MOHA) Methods and classes used to compute both molecular and spatial heterogeneity metrics
    // ===============================================================================================================================

	/**
	 * Default critical distance factor of 1.31 that is described in the MOHA tool publication. 
	 */
    public static double DEFAULT_CRIT_DISTANCE_FACTOR = 1.31;

    
	/**
	 * Default n-State threshold model of 3 representing a low, medium, and high categorical state for a biomarker. 
	 */
    public static int DEFAULT_N_STATE_THRESHOLD_MODEL = 3; 

	/**
	 * Default cell area of 443 square pixels and a cell radius of 11.875 pixels for a standard 20x magnification image.
	 */
	public static int DEFAULT_CELL_AREA = 443;

	/**
	 * The diverity types that are described in the MOHA Tool publication and are computed by this tool.
	 */
    public static enum Diversity_Type {Molecular, CellFamily, CellNeighbor, CellSocial, RandomizeOverLocations_CellFamily, CellCoordinationNumber, Molecular_Disparity};

    
	/**
	 * This is the maximum number of cell neighbors to report the cell coordination number frequency.
	 */
	protected int numCellNeighborsByFrequencyToReport = 10;

	/**
	 * Class cell data Table model.
	 */
	public class CellDataTableModel extends DefaultTableModel {
		
		protected static final long serialVersionUID = 1L;

		
		public String DEFAULT_BIOMARKER_COLNAME_TAG = "_Cell_Median";
	    
		public int MAX_NUM_CELLS = 2000000; // 2 million cells - required to prevent potential memory overflow
		
		public String colName_Sample_ID = "SAMPLE_ID";
	    public String colName_Slide_ID = "Slide_ID";
	    public String colName_Position_ID = "Position_ID";
	    public String colName_Cell_ID = "Cell_ID";
		public String colName_Cell_Center_X = "Cell_Center_X";  // Required for spatial metrics
		public String colName_Cell_Center_Y = "Cell_Center_Y";  // Required for spatial metrics
		public String colName_Cell_Area = "Cell_Area";  // Required for spatial metrics - alternative is to provide Cell_Radius
		public String colName_Cell_Radius = "Cell_Radius";  // Required for spatial metrics - alternative is to provide Cell_Area
		public String colName_Marker_States = "Marker_States";  // Required
	    
	    public int colIndx_Sample_ID = -1;
	    public int colIndx_Slide_ID = -1;
	    public int colIndx_Position_ID = -1;
	    public int colIndx_Cell_ID = -1;
		public int colIndx_Cell_Center_X = -1;
		public int colIndx_Cell_Center_Y = -1;
		public int colIndx_Cell_Radius = -1;
		public int colIndx_Marker_States = -1;
		
		public CellDataTableModel() {
			super();
		}
		
    	public List<String> getRequiredColNameList() {
        	List<String> requiredColNameList = new ArrayList<String>();
        	requiredColNameList.add(colName_Marker_States);
    		return requiredColNameList;
    	}
    	
    	public List<String> getAltColNameList() {
        	List<String> altColNameList = new ArrayList<String>();
        	altColNameList.add(colName_Sample_ID);
        	altColNameList.add(colName_Slide_ID);
        	altColNameList.add(colName_Position_ID);
        	altColNameList.add(colName_Cell_ID);
        	altColNameList.add(colName_Cell_Center_X);
        	altColNameList.add(colName_Cell_Center_Y);
        	altColNameList.add(colName_Cell_Area);
        	altColNameList.add(colName_Cell_Radius);
    		return altColNameList;
    	}
    	
    	public void parseAndValidateData(String cmsfileNameAndPath) throws Exception {

			if (hasRequiredSpatialColumns()) {
		    	int csi = -1;
		    	try {
			    	for(csi = 0;csi<getRowCount();csi++) {
			    		double x = Double.parseDouble((String)getValueAt(csi, findColumn(colName_Cell_Center_X)));
			    		setValueAt(x, csi, findColumn(colName_Cell_Center_X));
			    		double y = Double.parseDouble((String)getValueAt(csi, findColumn(colName_Cell_Center_Y)));
			    		setValueAt(y, csi, findColumn(colName_Cell_Center_Y));
			    		if (findColumn(colName_Cell_Area) >= 0) {
				    		double area = Double.parseDouble((String)getValueAt(csi, findColumn(colName_Cell_Area)));
				    		setValueAt(area, csi, findColumn(colName_Cell_Area));
			    		}
			    		if (findColumn(colName_Cell_Radius) >= 0) {
				    		double area = Double.parseDouble((String)getValueAt(csi, findColumn(colName_Cell_Radius)));
				    		setValueAt(area, csi, findColumn(colName_Cell_Radius));
			    		}
			    	}
		    	} catch (Exception e) {
		    		throw new Exception("Problems parsing data in row " + (csi + 1) + " in file " + cmsfileNameAndPath);
		    	}
				if (findColumn(colName_Cell_Radius) < 0) {
					//need to create Cell_Radius column from Cell_Area column data
			    	try {
				    	Double [] radiusData = new Double[getRowCount()];
				    	for(csi = 0;csi<getRowCount();csi++) {
				    		double area = (Double)getValueAt(csi, findColumn(colName_Cell_Area));
				    		radiusData[csi] = Math.sqrt(area / Math.PI);
				    	}
				    	addColumn(colName_Cell_Radius, radiusData);
			    	} catch (Exception e) {
			    		throw new Exception("Problem computing cell radius from Cell Area in row " + (csi + 1) + " in file " + cmsfileNameAndPath);
			    	}
				}
			}
			
    	}
    	
		public void updateColumnIndexes() {
			
			
			colIndx_Sample_ID = findColumn(colName_Sample_ID);
			colIndx_Slide_ID = findColumn(colName_Slide_ID);
			colIndx_Position_ID = findColumn(colName_Position_ID);
			colIndx_Cell_ID = findColumn(colName_Cell_ID);
			colIndx_Cell_Center_X = findColumn(colName_Cell_Center_X);
			colIndx_Cell_Center_Y = findColumn(colName_Cell_Center_Y);
			colIndx_Cell_Radius = findColumn(colName_Cell_Radius);
			colIndx_Marker_States = findColumn(colName_Marker_States);
		}
		
		public String getCellMarkerStates(int i) {
			return (String) getValueAt(i, colIndx_Marker_States);
		}
		
		public double getCellCenterX(int i) {
			return (Double) getValueAt(i, colIndx_Cell_Center_X);
		}
		
		public double getCellCenterY(int i) {
			return (Double) getValueAt(i, colIndx_Cell_Center_Y);
		}
		
		public double getCell_Radius(int i) {
			return (Double) getValueAt(i, colIndx_Cell_Radius);
		}
		

		public boolean hasRequiredSpatialColumns() {
			boolean requiredSpatialInfo = true;
			if (findColumn(colName_Cell_Center_X) < 0 || findColumn(colName_Cell_Center_Y) < 0)
				requiredSpatialInfo = false;
			if (findColumn(colName_Cell_Area) < 0 && findColumn(colName_Cell_Radius) < 0)
				requiredSpatialInfo = false;
			return requiredSpatialInfo;
		}
		
		public String getSampleID(int i) {
			if (colIndx_Sample_ID < 0) {
				return getSlideID(i) + "_"  + getPositionID(i);
			} else
				return (String) getValueAt(i, colIndx_Sample_ID);
		}
		
		public String getSlideID(int i) {
			if (colIndx_Slide_ID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_Slide_ID);
		}
		
		public String getPositionID(int i) {
			if (colIndx_Position_ID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_Position_ID);
		}
		
		public String getCellID(int i) {
			if (colIndx_Cell_ID < 0)
				return "";
			else
				return (String) getValueAt(i, colIndx_Cell_ID);
		}
	}
    


	/**
	 * Class for a state
	 */
	public class StateRec implements Comparable<StateRec> {
		
		public String stateID = "";
		public int rank = 0;
		public double numberOf = 0;
		public double frequency = 0;
		
		public StateRec() {

		}
		
		public StateRec(String stateID, int number) {
			this.stateID = stateID;
			this.numberOf = number;
		}

		public StateRec(int rank, int number) {
			this.rank = rank;
			this.numberOf = number;
		}
		

		public StateRec(StateRec o) {
			this.stateID = o.stateID;
			this.rank = o.rank;
			this.numberOf = o.numberOf;
		}
		
		public int compareTo(StateRec o) {
			if (o.numberOf < numberOf)
				return -1;
			else if (o.numberOf > numberOf)
				return 1;
			else
				return 0;
		}

	}
	

	/**
	 * Class for a diversity type.
	 */
	public class DiversityRec {
		
		public String diversityType = "";
		public double entropy = 0;
		public double numEntities = 0;
		public double heterogeneity = 0;
		public int numStates = 0;
		public int maxNumStates;
		public double avgStateIndex = 0;
		public int maxStateIndex = 0;

		public DiversityRec(String diversityType) {
			this.diversityType = diversityType;
		}
		
		public DiversityRec(Diversity_Type diversityType) {
			this(diversityType.toString());
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + ((diversityType == null) ? 0 : diversityType.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			DiversityRec other = (DiversityRec) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (diversityType == null) {
				if (other.diversityType != null)
					return false;
			} else if (!diversityType.equals(other.diversityType))
				return false;
			return true;
		}

		private MOHAtool getOuterType() {
			return MOHAtool.this;
		}

	}


	/**
	 * Returns a DiversityRec object with the computed entropy and heterogeneity. 
	 * The diversityType argument is used to identify the type of diversity that is being computed.
	 * The states argument is the distribution of entities, species, or things
	 * The maxNumStates argument is the maximum number of possible states that can exist in the system.
	 * <p>
	 * This method computes entropy of the states based on the Shannon-Weaver Diversity Index.
	 * The heterogeneity is computed by dividing the entropy by the log of the maximum number of states in the system.
	 *
	 * @param  diversityType  the diversity type or flavor
	 * @param  states         the distribution of the states
	 * @param  maxNumStates   the maximum number of states
	 * @return                the computed diversity metrics
	 * @see                   DiversityRec
	 * @see                   StateRec
	 * @see                   Diversity_Type
	 */
	public DiversityRec calcDiversity(String diversityType, Collection<StateRec> states, int maxNumStates) {

		DiversityRec diversity = new DiversityRec(diversityType);
		diversity.maxNumStates = maxNumStates;
		if (states == null || states.size() == 0)
			return diversity;

		diversity.numEntities = 0;
		diversity.numStates = 0;
		double sumStateIndexNum = 0;
		double sumNum = 0;
		for (StateRec state : states) {
			if (state.numberOf > 0) {
				diversity.numEntities += state.numberOf;
				diversity.numStates++;
				if (diversity.maxStateIndex < state.rank) {
					diversity.maxStateIndex = state.rank;
				}
				sumStateIndexNum += state.rank * state.numberOf;
				sumNum += state.numberOf;
			}
		}
		if (sumNum > 0) {
			diversity.avgStateIndex = sumStateIndexNum / sumNum;
		}

		// Compute the entropy of the system using the Shannon-Weaver Diversity Index
		diversity.entropy = 0;
		if (diversity.numEntities > 0) {
			for (StateRec state : states) {
				if (state.numberOf > 0) {
					double propValue = state.numberOf / (double) diversity.numEntities;
					if (propValue > 1e-30) {
						diversity.entropy += propValue * Math.log(propValue);
					}
				}
			}
		}
		diversity.entropy = Math.abs(diversity.entropy);
		
		// Compute the heterogeneity of the system. This value can range from zero to unity.
		if (diversity.maxNumStates > 1) {
			diversity.heterogeneity = diversity.entropy / Math.log(diversity.maxNumStates);
		} else {
			diversity.heterogeneity = 0;
		}

		return diversity;
	}



	/**
	 * Computes the molecular and spatial diversity metrics as described in the MOHA Tool publication.
	 * @param cellData                Inputted cell data table model.
	 * @param maxNumMolecularStates   The maximum number of molecular states for the system.
	 * @param criticalDistanceFactor  The critical distance factor to determine if two cells are neighbors.
	 * @return DiversityRec           Returns list of computed Diversity objects.
	 * @see CellDataTableModel
	 * @see Diversity_Type
	 */
	public List<DiversityRec> computeMOHAmetrics(CellDataTableModel cellData, int maxNumMolecularStates, double criticalDistanceFactor) throws Exception {

		cellData.updateColumnIndexes();
		
		int numCells = cellData.getRowCount();
		
		if (numCells == 0) {
			throw new Exception("Problem - No Cells (rows) in cell states data");
		}

		
		Map<String, StateRec> cellStateMap = new HashMap<String, StateRec>();
		for (int i = 0; i < cellData.getRowCount(); i++) {
			String cellMarkerState = cellData.getCellMarkerStates(i);
			StateRec state = cellStateMap.get(cellMarkerState);
			if (state == null) {
				cellStateMap.put(cellMarkerState, new StateRec(cellMarkerState, 1));
			} else {
				state.numberOf++;
			}
		}
		List<StateRec> cellStates = new ArrayList<StateRec>(cellStateMap.values());
		
		Collections.sort(cellStates);

    	for (int i=0;i<cellStates.size();i++) {
    		StateRec state = cellStates.get(i);
    		state.rank = i + 1;
    		state.frequency = state.numberOf / (double) numCells;
			// System.out.println(state.stateID + "\t" + state.numberOf + "\t" + state.frequency);
        }


		List<DiversityRec> mohaDiversityRecs = new ArrayList<DiversityRec>();
		
		int maxNumPossibleMolecularStates = Math.min(maxNumMolecularStates, numCells);
        DiversityRec molecularDiversity = calcDiversity(Diversity_Type.Molecular.toString(), cellStates, maxNumPossibleMolecularStates);
        mohaDiversityRecs.add(molecularDiversity);
        
    	DiversityRec molecularDisparityDiversity = new DiversityRec(Diversity_Type.Molecular_Disparity);
    	molecularDisparityDiversity.avgStateIndex = calcMolecularDisparity(cellStates);
        mohaDiversityRecs.add(molecularDisparityDiversity);
    	
        if (cellData.hasRequiredSpatialColumns()) {

	        int Nm = cellStates.size();
	        int Nc = 0;
	        for (int i = 0; i < Nm; i++) {
	        	Nc += cellStates.get(i).numberOf;
	        }
	        double [] Pmi = new double[Nm];
	        for (int i = 0; i < Nm; i++) {
	        	Pmi[i] = cellStates.get(i).numberOf / (double)Nc;
	        }
	
	        double sqrCritDistFactor = criticalDistanceFactor * criticalDistanceFactor;
	        
	        double [] fc = new double[numCellNeighborsByFrequencyToReport];
	        List<DiversityRec> diversityRecs = calcCellSpatialDiversity(cellData, sqrCritDistFactor, Nm,  Pmi, Nc, fc);
	        mohaDiversityRecs.addAll(diversityRecs);

	
	        DiversityRec cellSocialDiversity = calcCellSocialSpatialDiversity(cellStates, cellData, sqrCritDistFactor);
	        mohaDiversityRecs.add(cellSocialDiversity);

        }
        
         return mohaDiversityRecs;
	}
 
    

	/**
	 * Computes the molecular disparity for the inputted cell states.
	 * @param cellStates         Inputted cell molecular state frequency distribution.
	 * @see StateRec
	 */
	public double calcMolecularDisparity(List<StateRec> cellStates) {

		// DOUBLE precision does matter for this routine
		double molecularDisparity = 0;
		for (int i = 0; i < cellStates.size() - 1; i++) {
			StateRec pathwayStateRec_i = cellStates.get(i);
			char[] ciStateValues = pathwayStateRec_i.stateID.toCharArray();

			for (int j = i + 1; j < cellStates.size(); j++) {
				StateRec pathwayStateRec_j = cellStates.get(j);
				char[] cjStateValues = pathwayStateRec_j.stateID.toCharArray();

				int sqrStateDifference_i_and_j = 0;
				for (int k = 0; k < ciStateValues.length; k++) {
					int delta = cjStateValues[k] - ciStateValues[k];
					sqrStateDifference_i_and_j += delta * delta;
				}

				molecularDisparity += cellStates.get(i).frequency * cellStates.get(j).frequency * sqrStateDifference_i_and_j;
			}
		}

		return molecularDisparity;
	}
	
	
	/**
	 * Determines if two cells are considered neighbors.
	 * @param cellData            Inputted cell data table model.
	 * @param cell_i              Index of first cell i in cell data.
	 * @param cell_j              Index of second cell j in cell data.
	 * @param sqrCritDistFactor   The square of the critical distance factor to determine if two cells are neighbors.
	 * @return boolean            Returns true if the cell i and cell j are neighbors.
	 * @see CellDataTableModel
	 */
	public boolean isSpatialCellNeighbors(CellDataTableModel cellData, int cell_i, int cell_j, double sqrCritDistFactor) {

		if (!cellData.getPositionID(cell_i).equals(cellData.getPositionID(cell_j))  || !cellData.getSlideID(cell_i).equals(cellData.getSlideID(cell_j))) {
			return false;
		}
		
		double dx = cellData.getCellCenterX(cell_i) - cellData.getCellCenterX(cell_j);
		double dy = cellData.getCellCenterY(cell_i) - cellData.getCellCenterY(cell_j);
		//add the radius of each cell
		double dr = cellData.getCell_Radius(cell_i) + cellData.getCell_Radius(cell_j);
		
		if (dr > 0) {
			double normDist = (dx * dx + dy * dy) / (dr * dr);
			if (normDist <= sqrCritDistFactor)
				return true;
		}

		return false;
	}


	/**
	 * Calculates the spatial diversity for inputted cells including the cell family, cell neighbor, and cell social Diversity_Type
	 * @param cellData             Inputted cell data table model.
	 * @param sqrCritDistFactor    The square of the critical distance factor to determine if two cells are neighbors.
	 * @param Nm                   Number of molecular states.
	 * @param Pmi                  Frequency of molecular states.
	 * @param Nc                   Number of cells.
	 * @param fc                   Output of cell coordination (z) number fequency ditribution.
	 * @return DiversityRec        Returns list of computed Diversity objects.
	 * @see CellDataTableModel
	 * @see Diversity_Type
	 */
    public List<DiversityRec> calcCellSpatialDiversity(CellDataTableModel cellData, double sqrCritDistFactor, int Nm, double [] Pmi, int Nc, double [] fc) {

        Map<Integer, StateRec> cellFamilyEntitiesMap = new HashMap<Integer, StateRec>();
        Map<Integer, StateRec> cellNeighborEntitiesMap = new HashMap<Integer, StateRec>();
        Map<Integer, StateRec> cellCoordinationNumbersMap = new HashMap<Integer, StateRec>();
        
        int maxNumCellNeighbors = 0;
        for (int cell_i = 0; cell_i < cellData.getRowCount(); cell_i++) {

            Map<String, StateRec> neighborCellStates = new HashMap<String, StateRec>();
            int numCellNeighbors = 0;
            for (int cell_j = 0; cell_j < cellData.getRowCount(); cell_j++) {
                if (isSpatialCellNeighbors(cellData, cell_i, cell_j, sqrCritDistFactor)) {
                	if (cell_i != cell_j) {
	                    String cell_j_MarkerState = cellData.getCellMarkerStates(cell_j);
	                    StateRec neighborCellState = neighborCellStates.get(cell_j_MarkerState);
	                    if (neighborCellState == null) {
	                        neighborCellState = new StateRec(cell_j_MarkerState, 1);
	                        neighborCellStates.put(cell_j_MarkerState, neighborCellState);
	                    } else {
	                        neighborCellState.numberOf++;
	                    }
	                    numCellNeighbors++;
                	}
                }
            }

            maxNumCellNeighbors = Math.max(maxNumCellNeighbors, numCellNeighbors);

            //Group cells by how many cell neighbors they have (i.e. cell coordination number distribution)
            StateRec state = cellCoordinationNumbersMap.get(numCellNeighbors);
            if (state == null) {
            	state = new StateRec(numCellNeighbors, 1);
                cellCoordinationNumbersMap.put(numCellNeighbors, state);
            } else {
            	state.numberOf++;
            }
            
            //Group cells by cell states around touching cell measure i
            for (StateRec neighborCellState : neighborCellStates.values()) {
            	state = cellNeighborEntitiesMap.get((int)neighborCellState.numberOf);
                if (state == null) {
                	state = new StateRec((int)neighborCellState.numberOf, 1);
                    cellNeighborEntitiesMap.put((int)neighborCellState.numberOf, state);
                } else {
                	state.numberOf++;
                }
            }

            //Group cells by cell family (same cell state as cell measure i
            String cell_i_MarkerState = cellData.getCellMarkerStates(cell_i);
            
            int numNeighborsWithSameCellStatePlusOne = 1;
            StateRec neighborCellState = neighborCellStates.get(cell_i_MarkerState);
            if (neighborCellState != null) {
                numNeighborsWithSameCellStatePlusOne = (int)neighborCellState.numberOf + 1;
            }
            state = cellFamilyEntitiesMap.get(numNeighborsWithSameCellStatePlusOne);
            if (state == null) {
            	state = new StateRec(numNeighborsWithSameCellStatePlusOne, 1);
                cellFamilyEntitiesMap.put(numNeighborsWithSameCellStatePlusOne, state);
            } else {
            	state.numberOf++;
            }

        }

        List<StateRec> cellCoordinationNumbers = new ArrayList<StateRec>(cellCoordinationNumbersMap.values());
        // add one because a cell could have no neighbors
        int maxNumPossibleCellNeighborStates = maxNumCellNeighbors + 1;
        DiversityRec cellCoordinationNumberDiversity = calcDiversity(Diversity_Type.CellCoordinationNumber.toString(), cellCoordinationNumbers, maxNumPossibleCellNeighborStates);
        for(int i=0;i<numCellNeighborsByFrequencyToReport;i++) {
        	fc[i] = 0;
        }
        double sumNumEntities = 0;
        for(StateRec zState : cellCoordinationNumbers) {
        	if (zState.rank < numCellNeighborsByFrequencyToReport - 1) {
        		fc[zState.rank] += zState.numberOf;
        	} else {
        		fc[numCellNeighborsByFrequencyToReport - 1] += zState.numberOf;
        	}
        	sumNumEntities += zState.numberOf;
        }
        for(int i=0;i<numCellNeighborsByFrequencyToReport;i++) {
        	fc[i] /= sumNumEntities;
        }

        DiversityRec randomizeOverLocationsFamilyDiversity = computeRandomizedFamilyDiveristy(Diversity_Type.RandomizeOverLocations_CellFamily.toString(), Nm,  Pmi, Nc, cellCoordinationNumbers);
        
        List<StateRec> cellNeighborEntities = new ArrayList<StateRec>(cellNeighborEntitiesMap.values());
        DiversityRec cellNeighborDiversity = calcDiversity(Diversity_Type.CellNeighbor.toString(), cellNeighborEntities, maxNumPossibleCellNeighborStates);

        List<StateRec> cellFamilyEntities = new ArrayList<StateRec>(cellFamilyEntitiesMap.values());
        DiversityRec cellFamilyDiversity = calcDiversity(Diversity_Type.CellFamily.toString(), cellFamilyEntities, maxNumPossibleCellNeighborStates);

        List<DiversityRec> diversityRecs = new ArrayList<DiversityRec>();
        diversityRecs.add(cellNeighborDiversity);
        diversityRecs.add(cellFamilyDiversity);
        diversityRecs.add(cellCoordinationNumberDiversity);
        diversityRecs.add(randomizeOverLocationsFamilyDiversity);

        return diversityRecs;

    }    
  

	/**
	 * Calculates the cell social spatial diversity for the inputted cells.
	 * @param cellStates         Inputted cell molecular state frequency distribution.
	 * @param cellData           Inputted cell data table model.
	 * @param sqrCritDistFactor  The square of the critical distance factor to determine if two cells are neighbors.
	 * @return DiversityRec
	 * @see CellDataTableModel
	 * @see StateRec
	 * @see Diversity_Type
	 * @see DiversityRec
	 */
	public DiversityRec calcCellSocialSpatialDiversity(List<StateRec> cellStates, CellDataTableModel cellData, double sqrCritDistFactor) {

		Map<Integer, StateRec> cellSocialEntities = new HashMap<Integer, StateRec>();
		int numCells = 0;
		for (int i = 0; i < cellStates.size(); i++) {
			String selectedMarkerState = cellStates.get(i).stateID;

			List<Integer> selectedCells = new ArrayList<Integer>();
			for (int csi = 0; csi < cellData.getRowCount(); csi++) {
				if (selectedMarkerState.equals(cellData.getCellMarkerStates(csi))) {
					selectedCells.add(csi);
				}
			}

			class CellSocialLink implements Comparable<CellSocialLink> {
				int cellIndx;
				int groupIndx = -1;

				public CellSocialLink(int cellIndx, int groupIndx) {
					this.cellIndx = cellIndx;
					this.groupIndx = groupIndx;
				}

				public int compareTo(CellSocialLink o) {
					if (this.groupIndx > o.groupIndx)
						return 1;
					else if (this.groupIndx < o.groupIndx)
						return -1;
					else
						return 0;
				}
			}

			List<CellSocialLink> cellSocialLinks = new ArrayList<CellSocialLink>();
			for (int ci = 0; ci < selectedCells.size(); ci++) {
				cellSocialLinks.add(new CellSocialLink(ci, -1));
			}

			int numGroups = 0;
			for (int cgi = 0; cgi < cellSocialLinks.size(); cgi++) {
				CellSocialLink cellGroupLink_i = cellSocialLinks.get(cgi);
				int cell_i = selectedCells.get(cellGroupLink_i.cellIndx);

				for (int cgj = 0; cgj < cellSocialLinks.size(); cgj++) {

					if (cgj == cgi)
						continue;

					CellSocialLink cellGroupLink_j = cellSocialLinks.get(cgj);

					if (cellGroupLink_i.groupIndx != -1 && cellGroupLink_i.groupIndx == cellGroupLink_j.groupIndx) {
						continue;
					}

					int cell_j = selectedCells.get(cellGroupLink_j.cellIndx);

					if (isSpatialCellNeighbors(cellData, cell_i, cell_j, sqrCritDistFactor)) {

						if (cellGroupLink_i.groupIndx == -1 && cellGroupLink_j.groupIndx == -1) {
							cellGroupLink_i.groupIndx = numGroups;
							cellGroupLink_j.groupIndx = numGroups;
							numGroups++;
						} else if (cellGroupLink_i.groupIndx == -1) {
							cellGroupLink_i.groupIndx = cellGroupLink_j.groupIndx;
						} else if (cellGroupLink_j.groupIndx == -1) {
							cellGroupLink_j.groupIndx = cellGroupLink_i.groupIndx;

						} else if (cellGroupLink_i.groupIndx < cellGroupLink_j.groupIndx) {
							int deleteGroupIndx = cellGroupLink_j.groupIndx;
							int combineGroupIndx = cellGroupLink_i.groupIndx;
							for (int ci = 0; ci < cellSocialLinks.size(); ci++) {
								if (cellSocialLinks.get(ci).groupIndx == deleteGroupIndx)
									cellSocialLinks.get(ci).groupIndx = combineGroupIndx;
							}
							for (int gi = deleteGroupIndx + 1; gi < numGroups; gi++) {
								for (int ci = 0; ci < cellSocialLinks.size(); ci++) {
									if (cellSocialLinks.get(ci).groupIndx == gi)
										cellSocialLinks.get(ci).groupIndx--;
								}
							}
							numGroups--;
						} else if (cellGroupLink_i.groupIndx > cellGroupLink_j.groupIndx) {
							int deleteGroupIndx = cellGroupLink_i.groupIndx;
							int combineGroupIndx = cellGroupLink_j.groupIndx;
							for (int ci = 0; ci < cellSocialLinks.size(); ci++) {
								if (cellSocialLinks.get(ci).groupIndx == deleteGroupIndx)
									cellSocialLinks.get(ci).groupIndx = combineGroupIndx;
							}
							for (int gi = deleteGroupIndx + 1; gi < numGroups; gi++) {
								for (int ci = 0; ci < cellSocialLinks.size(); ci++) {
									if (cellSocialLinks.get(ci).groupIndx == gi)
										cellSocialLinks.get(ci).groupIndx--;
								}
							}
							numGroups--;
						} else if (cellGroupLink_i.groupIndx == cellGroupLink_j.groupIndx) {
							System.out.println("SHOULD NEVER REACH THIS POINT");
						}

					}

				}
			}


			// Assign single cells into one group
			List<String> singleCellGroups = new ArrayList<String>();
			for (int ci = 0; ci < cellSocialLinks.size(); ci++) {
				CellSocialLink cellGroupLink = cellSocialLinks.get(ci);
				if (cellGroupLink.groupIndx == -1) {
					String cell_index = selectedCells.get(cellGroupLink.cellIndx).toString();
					singleCellGroups.add(cell_index);
				}
			}

			Map<Integer, StateRec> stateMap = new HashMap<Integer, StateRec>();
			if (singleCellGroups.size() > 0) {
				StateRec state = new StateRec();
				state.rank = 1;
				state.numberOf = singleCellGroups.size();
				stateMap.put(state.rank, state);
			}
			for (int groupIndx = 0; groupIndx < numGroups; groupIndx++) {
				int numCellsInCluster = 0;
				for (int ci = 0; ci < cellSocialLinks.size(); ci++) {
					CellSocialLink cellGroupLink = cellSocialLinks.get(ci);
					if (cellGroupLink.groupIndx == groupIndx) {
						numCellsInCluster++;
					}
				}

				StateRec state = stateMap.get(numCellsInCluster);
				if (state != null) {
					state.numberOf++;
				} else {
					state = new StateRec();
					state.rank = numCellsInCluster;
					state.numberOf = 1;
					stateMap.put(state.rank, state);
				}

			}

			List<StateRec> cellSocialStates = new ArrayList<StateRec>(stateMap.values());

			numCells += selectedCells.size();
			for (StateRec cellSocialState : cellSocialStates) {
				StateRec mapClusterRec = cellSocialEntities.get(cellSocialState.rank);
				if (mapClusterRec == null) {
					cellSocialEntities.put(cellSocialState.rank, new StateRec(cellSocialState));
				} else {
					mapClusterRec.numberOf += cellSocialState.numberOf;
				}
			}

		}

        
        int maxPossibleStates_All = 1;
        if (numCells > 1) {
            maxPossibleStates_All =  (int)((Math.sqrt(8 * numCells + 1) - 1) / 2);
        }
		DiversityRec cellSocialDiversity = calcDiversity(Diversity_Type.CellSocial.toString(), cellSocialEntities.values(), maxPossibleStates_All);

		return cellSocialDiversity;
	}
    
	/**
	 * Computes a Randomized Family Diveristy metric for a regular lattice of cells with a coordination number z.
	 * @param diversityType Diversity_Type string to assign to the computed diversity
	 * @param Nm            Number of molecular states.
	 * @param Pmi           frequency of molecular states
	 * @param Nc            number of cells
	 * @param z             cell coordination number assuming a regular lattice
	 * @return DiversityRec
	 * @see Diversity_Type
	 * @see DiversityRec
	 */
	public DiversityRec computeRandomizedFamilyDiveristy(String diversityType, int Nm, double[] Pmi, int Nc, int z) {

		List<StateRec> zDist = new ArrayList<StateRec>();
		StateRec clusterRec = new StateRec();
		clusterRec.rank = z;  // regular lattice coordination number
		clusterRec.numberOf = 1;
		zDist.add(clusterRec);
		return computeRandomizedFamilyDiveristy(diversityType, Nm, Pmi, Nc, zDist);
	}

	/**
	 * Computes a Randomized Family Diveristy metric as described in the MOHA Tool publication. 
	 * @param diversityType Diversity_Type string to assign to the computed diversity
	 * @param Nm            Number of molecular states.
	 * @param Pmi           frequency of molecular states
	 * @param Nc            number of cells
	 * @param zStates       frequency distribution of cell coordination (z) numbers
	 * @return DiversityRec
	 * @see Diversity_Type
	 * @see DiversityRec
	 */
	public DiversityRec computeRandomizedFamilyDiveristy(String diversityType, int Nm, double[] Pmi, int Nc, List<StateRec> zStates) {

		int maxZ = 0;
		for (StateRec zState : zStates) {
			if (maxZ < zState.rank)
				maxZ = zState.rank;
		}

		double[] Pfk = new double[maxZ + 1];
		for (int k = 0; k <= maxZ; k++) {
			Pfk[k] = 0;
			for (StateRec zState : zStates) {
				if (zState.numberOf > 0) {
					int z = zState.rank;
					double Pfkz = 0;
					if (k == 0) {
						for (int i = 0; i < Nm; i++) {
							Pfkz += Pmi[i] * Math.pow((1.0 - Pmi[i]), z);
						}

					} else if (k == z) {
						for (int i = 0; i < Nm; i++) {
							Pfkz += Pmi[i] * Math.pow(Pmi[i], k);
						}
					} else {
						for (int i = 0; i < Nm; i++) {
							Pfkz += Pmi[i] * Math.pow((1.0 - Pmi[i]), (z - k)) * Math.pow(Pmi[i], k);
						}
						double numPerm = 1;
						int gz = z;
						int gk = k;
						while (gz > z - k | gk > 0) {
							if (gz > z - k) {
								numPerm *= gz;
								gz--;
							}
							if (gk > 0) {
								numPerm /= gk;
								gk--;
							}
						}
						Pfkz *= numPerm;

					}
					Pfk[k] += zState.numberOf * Pfkz;
				}
			}

		}

		List<StateRec> states = new ArrayList<StateRec>();
		for (int k = 0; k <= maxZ; k++) {
			StateRec state = new StateRec();
			state.rank = k;
			state.numberOf = Pfk[k];
			states.add(state);
		}

		int numPossibleStates = maxZ + 1; // add one because a cell could have no neighbors
		DiversityRec gridFamilyDiversity = calcDiversity(diversityType, states, numPossibleStates);

		return gridFamilyDiversity;
	}
	    
	/**
	 * Returns the diversity object from a list of diversity objects of the requested diversity type or null otherwise.
	 * @param diversityRecs  List of diversity objects.
	 * @param diversityType  Diversity_Type string. 
	 * @return DiversityRec
	 * @see Diversity_Type
	 * @see DiversityRec
	 */
	protected DiversityRec getDiversityType(List<DiversityRec> diversityRecs, String diversityType) {
		int indx = diversityRecs.indexOf(new DiversityRec(diversityType));
		if (indx < 0)
			return null;
		else
			return diversityRecs.get(indx);
	}


	/**
	 * Prints diversity metrics to an output file
	 * @param sampleID              The sample ID for the inputted list of diversity objects
	 * @param mohaContext           The context for the computed diversity (e.g. pathway, gene, or marker set used to define the cell state.
	 * @param diversityRecs         List of diversity objects to print out in a single row
	 * @param numCells              The number of cells used to compute the list of diversity objects
	 * @param maxNumCellStates      Maximum number of molecular cell states.
	 * @param df                    DecimalFormat used to format output of diversity metrics
	 * @param outputFileNameAndPath Output filename and path to print diversity metrics to.
	 * @param append                Set to true to append output to an existing file or false to overwrite file.
	 * @see DiversityRec
	 */
	protected void printDiversityMetrics(String sampleID, String mohaContext, List<DiversityRec> diversityRecs, int numCells, int maxNumCellStates, DecimalFormat df, String outputFileNameAndPath, boolean append) throws Exception {


		int maxNumPossibleMolecularStates = Math.min(maxNumCellStates, numCells);
		
		DiversityRec molecularDiversity = getDiversityType(diversityRecs, Diversity_Type.Molecular.toString());
		DiversityRec molecularDisparityDiversity = getDiversityType(diversityRecs, Diversity_Type.Molecular_Disparity.toString());
		
		DiversityRec cellNeighborDiversity = getDiversityType(diversityRecs, Diversity_Type.CellNeighbor.toString());
		DiversityRec cellFamilyDiversity = getDiversityType(diversityRecs, Diversity_Type.CellFamily.toString());
		DiversityRec cellCoordinationNumberDiversity = getDiversityType(diversityRecs, Diversity_Type.CellCoordinationNumber.toString());
		DiversityRec randomizeOverLocationFamilyDiversity = getDiversityType(diversityRecs, Diversity_Type.RandomizeOverLocations_CellFamily.toString());
		DiversityRec cellSocialDiversity = getDiversityType(diversityRecs, Diversity_Type.CellSocial.toString());
			

        boolean printHeader = !append;
        
      	File outputFile = new File(outputFileNameAndPath);
    	if (!outputFile.exists()) {
    		printHeader = true;
    	}	
    	
  
        //Note that new FileWriter Append flag is false
        PrintWriter pw = new PrintWriter(new FileWriter(outputFileNameAndPath, append));
        
        
        if (printHeader) {
	        String outputHeaderLine = "SAMPLE_ID"  + "\t" + "MOHA_CONTEXT";
	        outputHeaderLine += "\t" + "NumCells";
	        outputHeaderLine += "\t" + "AvgCellCoordinationNumber" + "\t" + "MaxCellCoordinationNumber";
	        outputHeaderLine += "\t" + "CellCoordinationNumber_Entropy" + "\t" + "CellCoordinationNumber_Heterogeneity";
	        outputHeaderLine += "\t" + "MaxNumStates" + "\t" + "NumObsStates" + "\t" + "Molecular_Disparity";
	        outputHeaderLine += "\t" + "Molecular_Entropy" + "\t" + "Molecular_Heterogeneity";
	        outputHeaderLine += "\t" + "CellFamily_Entropy" + "\t" + "CellFamily_Heterogeneity";
	        outputHeaderLine += "\t" + "CellNeighbor_Entropy" + "\t" + "CellNeighbor_Heterogeneity";
	        outputHeaderLine += "\t" + "CellSocial_Entropy" + "\t" + "CellSocial_Heterogeneity";
	        outputHeaderLine += "\t" + "RandomizeOverLocations_CellFamily_Entropy" + "\t" + "RandomizeOverLocations_CellFamily_Heterogeneity";


        	pw.println(outputHeaderLine);
        }

		String dataLineSampleDiversity = "";

		dataLineSampleDiversity += sampleID + "\t" + mohaContext;
		dataLineSampleDiversity += "\t" + numCells;
		if (cellCoordinationNumberDiversity != null) {
			dataLineSampleDiversity += "\t" + df.format(cellCoordinationNumberDiversity.avgStateIndex) + "\t" + cellCoordinationNumberDiversity.maxStateIndex;
			dataLineSampleDiversity += "\t" + df.format(cellCoordinationNumberDiversity.entropy) + "\t" + df.format(cellCoordinationNumberDiversity.heterogeneity);
		} else {
			dataLineSampleDiversity += "\t" + "NA" + "\t" + "NA";
			dataLineSampleDiversity += "\t" + "NA" + "\t" + "NA";
		}

		dataLineSampleDiversity += "\t" + maxNumPossibleMolecularStates + "\t" + molecularDiversity.numStates + "\t" + df.format(molecularDisparityDiversity.avgStateIndex);
		dataLineSampleDiversity += "\t" + df.format(molecularDiversity.entropy) + "\t" + df.format(molecularDiversity.heterogeneity);

		if (cellFamilyDiversity != null) {
			dataLineSampleDiversity += "\t" + df.format(cellFamilyDiversity.entropy) + "\t" + df.format(cellFamilyDiversity.heterogeneity);
		} else {
			dataLineSampleDiversity += "\t" + "NA" + "\t" + "NA";
		}
		if (cellNeighborDiversity != null) {
			dataLineSampleDiversity += "\t" + df.format(cellNeighborDiversity.entropy) + "\t" + df.format(cellNeighborDiversity.heterogeneity);
		} else {
			dataLineSampleDiversity += "\t" + "NA" + "\t" + "NA";
		}
		if (cellSocialDiversity != null) {
			dataLineSampleDiversity += "\t" + df.format(cellSocialDiversity.entropy) + "\t" + df.format(cellSocialDiversity.heterogeneity);
		} else {
			dataLineSampleDiversity += "\t" + "NA" + "\t" + "NA";
		}
		if (randomizeOverLocationFamilyDiversity != null) {
			dataLineSampleDiversity += "\t" + df.format(randomizeOverLocationFamilyDiversity.entropy) + "\t" + df.format(randomizeOverLocationFamilyDiversity.heterogeneity);
		} else {
			dataLineSampleDiversity += "\t" + "NA" + "\t" + "NA";
		}


		pw.println(dataLineSampleDiversity);
		pw.flush();
		pw.close();
	}

}






