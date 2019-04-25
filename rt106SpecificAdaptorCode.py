# Copyright (c) General Electric Company, 2019.  All rights reserved.

# Rt 106

# Multi-Omics Heterogeneity Analysis

import os, glob, uuid, time, json, string, logging, subprocess

# function: run_algorithm() -- Python function for marshalling your data and running your algorithm.
# parameters:
#   datastore: object to be used when interacting with the Data Store
#   context:  A JSON structure that should contain all the inputs and parameters your algorithm needs.
def run_algorithm(datastore, context):

    logging.info('run_algorithm: %r' % context)

    for f in glob.glob('/rt106/data/*'):
        os.remove(f)

    # Compute biomarker thresholds.
    cell_quant_path = datastore.get_pathology_result_image_path(context['slide'], context['region'], context['branch'],  'Quant')
    if (type(cell_quant_path) == "int" and cell_quant_path > 200):
        status = "ERROR_QUANT_FILE_NOT_FOUND"
        return { 'result' : {}, 'status' : cell_quant_path }

    quant_csv = 'quant_%s.csv' % context['region']
    instance_status = datastore.get_instance(cell_quant_path, '/rt106/data', quant_csv, 'csv')
    if (instance_status != 200):
        status = "ERROR_QUANT_FILE_NOT_FOUND"
        return { 'result' : {}, 'status' : status }

    MOHA_path = datastore.get_pathology_result_path(context['slide'], context['region'], context['branch'], 'MOHA')
    moha_thresholds_csv = 'quant_%s.csv.thresholds.txt' % context['region']
    moha_markerstates_csv = 'quant_%s.csv.MarkerStates.txt' % context['region']
    moha_markerindex_csv = 'quant_%s.csv.MarkerIndex.txt' % context['region']
    moha_heterogeneity_csv = 'out_moha_%s.txt' % context['region']
    #output_file_thresholds = '/rt106/input/%s' % moha_thresholds_csv

    # Code for calling algorithm.
    try:
        run_threshold_algorithm = 'java -jar MOHAtool.jar -computeThresholds=data/quant_%s.csv -biomarkerMetricTag=_Cell_Mean' % context['region']
        logging.info('run Algorithm: %r' % run_threshold_algorithm)
        subprocess.check_call(run_threshold_algorithm,shell=True)
    except subprocess.CalledProcessError, e:
        logging.error('%d - %s' % (e.returncode, e.cmd))
        status = "EXECUTION_FINISHED_ERROR"
        result_context = {}
        return { 'result' : result_context, 'status' : status }


    # Based on the success or failure of your algorithm execution, set the 'status' variable to an appropriate string.
    status = "EXECUTION_FINISHED_SUCCESS"

    # Store result files and obtain the URIs that refer to those files.
    response_threshold_upload = datastore.post_instance(MOHA_path, '/rt106/data', moha_thresholds_csv, 'csv', context['force'])

    if response_threshold_upload == 403:
        status = "EXECUTION_ERROR"

    if status == "EXECUTION_FINISHED_SUCCESS":
        # Calculate cell states.
        try:
            run_cell_states_algorithm = 'java -jar MOHAtool.jar -computeCellStates=data/quant_%s.csv -thresholdFile=data/quant_%s.csv.thresholds.txt' % (context['region'],context['region'])
            logging.info('run Algorithm: %r' % run_cell_states_algorithm)
            subprocess.check_call(run_cell_states_algorithm,shell=True)
        except subprocess.CalledProcessError, e:
            logging.error('%d - %s' % (e.returncode, e.cmd))
            status = "EXECUTION_FINISHED_ERROR"
            result_context = {}
            return { 'result' : result_context, 'status' : status }
        else:
            # Store result files and obtain the URIs that refer to those files.
            response_markerstates_upload = datastore.post_instance(MOHA_path, '/rt106/data', moha_markerstates_csv, 'csv', context['force'])
            response_markerindex_upload = datastore.post_instance(MOHA_path, '/rt106/data', moha_markerindex_csv, 'csv', context['force'])


    if status == "EXECUTION_FINISHED_SUCCESS":
        # Calculate heterogeneity metrics.
        try:
            run_heterogeneity_algorithm = 'java -jar MOHAtool.jar -computeHeterogeneity=data/quant_%s.csv.MarkerStates.txt -outputFile=data/out_moha_%s.txt -append=false' % (context['region'],context['region'])
                #'java -jar MOHAtool.jar -computeCellStates=example/quant_%s.csv -thresholdFile=example/quant_%s.csv.thresholds.txt' % (context['region'],context['region'])
            logging.info('run Algorithm: %r' % run_heterogeneity_algorithm)
            subprocess.check_call(run_heterogeneity_algorithm,shell=True)
        except subprocess.CalledProcessError, e:
            logging.error('%d - %s' % (e.returncode, e.cmd))
            status = "EXECUTION_FINISHED_ERROR"
            result_context = {}
            return { 'result' : result_context, 'status' : status }
        else:
            # Store result files and obtain the URIs that refer to those files.
            response_heterogeneity_upload = datastore.post_instance(MOHA_path, '/rt106/data', moha_heterogeneity_csv, 'csv', context['force'])

    # Parse heterogeneity metrics file to pull out key metrics.  Return the metrics and save them to the database.
    if status == "EXECUTION_FINISHED_SUCCESS":
        hetfile = open('/rt106/data/'+moha_heterogeneity_csv, "r")
        hetcolumns = hetfile.readline()
        hetvalues = hetfile.readline()
        hetfile.close()
        hetcol_list = hetcolumns.strip().split('\t')
        hetval_list = hetvalues.strip().split('\t')
        hetdict = {}
        for i in range(len(hetcol_list)):
            hetdict[hetcol_list[i]] = hetval_list[i]

    # 5.    Create JSON structure containing results.
    nuclear_image_path = datastore.get_pathology_primary_path(context['slide'], context['region'], 'DAPI')
    result_context = {
        "nuclearImage" : nuclear_image_path,
        "CellFamily_Heterogeneity" : hetdict['CellFamily_Heterogeneity'],
        "CellSocial_Heterogeneity" : hetdict['CellSocial_Heterogeneity'],
        "Molecular_Heterogeneity" : hetdict['Molecular_Heterogeneity']
    }

    return { 'result' : result_context, 'status' : status }
