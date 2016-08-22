""" 
Fusion PQ (Performance Quality) analysis program. The following
program takes in a modified combined LIS & PCR file and performs
checks on the data to see if it meets specifications set by 
a guideline
"""

import pandas as pd
import numpy as np
import os

class FusionPQ():

    def __init__(self, dframe, assay_type, pos_label, neg_label, *args, **kwargs):

            self.dframe = dframe
            self.fix_empty()

            if assay_type == "P 1/2/3/4":

                self.settings = {
                    "assay_type":assay_type,
                    "pos_ctrl":'103101.*',
                    "neg_ctrl":'101111.*',
                    "pos_label":pos_label,
                    "neg_label":neg_label,
                }

                self.thresholds = {
                    "POS":{
                        "FAM_RFU_THRESHOLD_MIN":1200,
                        "FAM_RFU_THRESHOLD_MAX":None,
                        "HEX_RFU_THRESHOLD_MIN":2000,
                        "HEX_RFU_THRESHOLD_MAX":None,
                        "ROX_RFU_THRESHOLD_MIN":1500,
                        "ROX_RFU_THRESHOLD_MAX":None,
                        "RED647_RFU_THRESHOLD_MIN":400,
                        "RED647_RFU_THRESHOLD_MAX":None,
                        "IC_RFU_THRESHOLD_MIN":None,
                        "IC_RFU_THRESHOLD_MAX":None                     
                    },
                    "NEG":{
                        "FAM_RFU_THRESHOLD_MIN":None,
                        "FAM_RFU_THRESHOLD_MAX":500,
                        "HEX_RFU_THRESHOLD_MIN":None,
                        "HEX_RFU_THRESHOLD_MAX":500,
                        "ROX_RFU_THRESHOLD_MIN":None,
                        "ROX_RFU_THRESHOLD_MAX":500,
                        "RED647_RFU_THRESHOLD_MIN":None,
                        "RED647_RFU_THRESHOLD_MAX":350,
                        "IC_RFU_THRESHOLD_MIN":2000,
                        "IC_RFU_THRESHOLD_MAX":None                     
                    }
                }


    def fix_empty(self):
         self.dframe['Valid/Invalid for IC'] = self.dframe['Valid/Invalid for IC'].fillna("Invalid")

    def run_pq(self):

        POS_CTRL = self.settings['pos_ctrl']
        NEG_CTRL = self.settings['neg_ctrl']
        POS_LBL = self.settings['pos_label']
        NEG_LBL = self.settings['neg_label']

        self.set_labels(POS_LBL, NEG_LBL, POS_CTRL, NEG_CTRL)
        self.check_false_calls()
        self.check_validity()
        self.check_pq_thresholds()
        self.overall_validity()

        self.dframe.to_excel(os.path.join(os.getcwd(), "tester.xlsx"))
        return self.dframe

    def set_labels(self, pos_label, neg_label, pos_ctrl, neg_ctrl):
        """ 
        Creates a temporary column that allows us to categorize the
        samples as a positive sample, a negative sample, or controls.

        Args:
            pos_label - a string pattern that is shared with all positive samples (str)
            neg_label - a string pattern that is shared with all negative samples (str)
            pos_ctrl - a regex pattern that corresponds to positive controls (str)
            neg_ctrl - a regex pattern that corresponds to negative controls (str)
        Returns:
            None
        """
        self.dframe['SAMPLE Category'] = self.dframe['Specimen Barcode'].apply(lambda row: self.__set_labels_helper(row, pos_label, neg_label, pos_ctrl, neg_ctrl))
        
    def __set_labels_helper(self, row, pos_label, neg_label, positive_ctrl, negative_ctrl):
        """ The column marker logic for set_labels method """

        import re

        # pos_pattern = re.compile(r'%s' % positive_ctrl)
        # neg_pattern = re.compile(r'%s' % negative_ctrl)

        if row.upper() in pos_label:
    
            return "POS"
        elif row.upper() in neg_label:
 
            return "NEG"
        else:
            return ""

        # if row.lower().find(pos_label) > -1 or pos_pattern.search(row.lower()):
        #     return "POS"
        # elif row.lower().find(neg_label) > -1 or neg_pattern.search(row.lower()):
        #     return "NEG"
        # else:
        #     return ""

    def check_validity(self):
        """ Checks the validity status of each channel.
        This only checks if the samples are considered 'valid' which is
        different from the samples passing PQ criterias. The validity status
        will be stored in a temporary column.

        Samples are classified as follows:
        (1) If the interpretation of a channel is invalid, sample is INVALID.
        (2) If a 'positive' sample is positive and IC is valid, sample is VALID.
        (3) If a 'positive' sample is negative and IC is valid, sample is VALID.
        (4) If a 'positive' sample is positive and IC is invalid, sample is VALID. *
        (5) If a 'positive' sample is negative and IC is invalid, sample is INVALID.
        (6) If a 'negative' sample is positive and IC is invalid, sample is INVALID.
        (7) If a 'negative' sample is positive and IC is valid, sample is VALID.
        (8) If a 'negative' sample is negative and IC is valid, sample is VALID.
        (9) If a 'negative' sample is negative and IC is invalid, sample is INVALID.

        * Not a typo, invalid IC's can be overshadowed by a high positive value.
        """

        columns_to_validate = {
            "FAM": "POS/NEG/Invalid for HPIV-1",
            "HEX": "POS/NEG/Invalid for HPIV-2",
            "RED647": "POS/NEG/Invalid for HPIV-4",
            "ROX": "POS/NEG/Invalid for HPIV-3"
        }

        for channel, column_name in columns_to_validate.items():
            results_column = "Validity for {0}".format(column_name)
            self.dframe[results_column] = self.dframe.apply(self.__check_validity_helper, axis=1, column_type=column_name)

    def __check_validity_helper(self, row, column_type):
        """ The column validity helper for check validity method 

        Args:
            row - the row to modify (pandas obj)
            column_type - the specific channel that's being interpreted. (str)
        Returns:
            'Valid' if the data is valid. (str)
            'Invalid' if the data is invalid. (str)
        """

        import re

        # Temp hotfix
        if row[column_type] == "":
            return "INVALID"

        channel_examine = row[column_type]

        # Very strange... sometime loads as float
        channel_examine = str(channel_examine)

        # Helps differentiate between positive or negative samples to apply validity logic.
        sample_category = row['SAMPLE Category']

        ic_column = row['Valid/Invalid for IC']

        label_for_validity = re.compile(r'\bvalid\b')

        # Sample type 'positive'
        if sample_category.lower().find("pos") > -1:
            if channel_examine.lower().find("pos") > -1:
                return "VALID"
            elif channel_examine.lower().find("neg") > -1 and label_for_validity.search(ic_column.lower()):
                return "VALID"
            else:
                return "INVALID"
        # Sample type 'negative'
        elif sample_category.lower().find("neg") > -1:
            if channel_examine.lower().find("neg") > -1 and label_for_validity.search(ic_column.lower()):
                return "VALID"
            elif channel_examine.lower().find("pos") > -1 and label_for_validity.search(ic_column.lower()):
                return "VALID"
            else:
                return "INVALID"
        else:
                return "INVALID"

        return "INVALID"

    def check_false_calls(self):
        """ 
        Has to be called after labels are set.

        Checks for False Positives or False Negatives. Creates a temporary
        column to hold the results for false positives or false negatives. If the
        data does not contain false pos/false neg, it will have no status marked in this column.

        Args:
            None
        Returns:
            None
        """

        columns_to_examine_false = {
            "FAM": "POS/NEG/Invalid for HPIV-1",
            "HEX": "POS/NEG/Invalid for HPIV-2",
            "RED647": "POS/NEG/Invalid for HPIV-4",
            "ROX": "POS/NEG/Invalid for HPIV-3"         
        }

        for channel, column_name in columns_to_examine_false.items():
            results_column = "TRUTHS table for {0}".format(column_name)
            self.dframe[results_column] = self.dframe.apply(self.__check_false_calls_helper, axis=1, column_type=column_name)

    def __check_false_calls_helper(self, row, column_type):
        """ Helper method for check false calls
        Args:
            row - the specific row to analyze (pandas obj)
            column_type - the specific channel that's being interpreted. (str)
        Returns:
            'False Negative' if result is false negative (str)
            'False Positive' if result is false positive (str)
        """

        sample_barcode = row['SAMPLE Category']
        channel_column = row[column_type]

        if sample_barcode.lower().find("POS") > -1:
            if channel_column.lower().find("neg") > -1:
                return "False Negative"
            return ""
        elif sample_barcode.lower().find("NEG") > -1:
            if channel_column.lower().find("pos") > -1:
                return "False Positive"
            return ""
        else:
            return ""

    def check_pq_thresholds(self):

        """ Checks if the RFU Ranges meet the required
        threshold set by the PQ guidelines. Creates a temporary
        column to hold the status of 'PASS' or 'FAIL'
        """

        columns_to_pq = {
            "FAM": "FAM Rounded RFU Range (HPIV-1)",
            "HEX": "HEX Rounded RFU Range (HPIV-2)",
            "RED647": "RED647 Rounded RFU Range (HPIV-4)",
            "ROX": "ROX Rounded RFU Range (HPIV-3)",
            "IC":"IC Rounded RFU Range"
        }

        for chan, channel_col in columns_to_pq.items():
            results_column = "PQ Threshold for {0}".format(channel_col)
            self.dframe[results_column] = self.dframe.apply(self.__pq_threshold_helper, channel_type=chan, channel_col_name=channel_col, axis=1)

    def __pq_threshold_helper(self, row, channel_type, channel_col_name):   
        """ Helper method for check_pq_threshold
        Args:
            row - the row to be analyzed (pd obj)
            channel_type - the channel being analyzed (str)
            channel_col_name - the actual column with RFU data (str)
        Returns:
            "PASS" if PQ is met (str)
            "FAIL" if PQ is not met (str)
        """

        sample_category = row['SAMPLE Category']

        thresh_cat = {}

        if sample_category.lower().find("pos") > -1:
            thresh_cat = self.thresholds['POS']
        else:
            thresh_cat = self.thresholds['NEG']

        THRESHOLD_SETTINGS = {'FAM':[thresh_cat["FAM_RFU_THRESHOLD_MIN"], thresh_cat["FAM_RFU_THRESHOLD_MAX"]],'HEX': [thresh_cat["HEX_RFU_THRESHOLD_MIN"], thresh_cat["HEX_RFU_THRESHOLD_MAX"]],'ROX':[thresh_cat["ROX_RFU_THRESHOLD_MIN"], thresh_cat["ROX_RFU_THRESHOLD_MAX"]],'RED647':[thresh_cat["RED647_RFU_THRESHOLD_MIN"], thresh_cat["RED647_RFU_THRESHOLD_MAX"]],'IC': [thresh_cat["IC_RFU_THRESHOLD_MIN"], thresh_cat["IC_RFU_THRESHOLD_MAX"]]}

        # Determine the comparison cases

        # Comparison Modes
        # 0 = checks if RFU range is greater than min
        # 1 = checks if RFU range is less than max
        # 2 = checks if RFU range is between min and max
        comparison_mode = 2
        channel_settings = THRESHOLD_SETTINGS[channel_type]
        
        if channel_settings[0] == None and channel_settings[1] == None:
                return "PASS"
        
        for number, check_none in enumerate(channel_settings):
                if check_none is not None:
                        comparison_mode = number

        # Quick hacky fix for postive samples with no data
        rfu_range = row[channel_col_name]
        
        if rfu_range == "-":
                rfu_range = 0
        else:
                rfu_range = float(rfu_range)

        if comparison_mode == 0:
                if rfu_range > channel_settings[0]:
                        return "PASS"
        elif comparison_mode == 1:
                if rfu_range < channel_settings[1]:
                        return "PASS"
        else:
                if channel_settings[0] < rfu_range < channel_settings[1]:
                        return "PASS"
                
        return "FAIL"

    def overall_validity(self):
        self.dframe['PQ RESULTS'] = self.dframe.apply(self.__overall_validity_helper, axis=1)

    def __overall_validity_helper(self, row):

        columns_to_validate = {
            "FAM": "Validity for POS/NEG/Invalid for HPIV-1",
            "HEX": "Validity for POS/NEG/Invalid for HPIV-2",
            "RED647": "Validity for POS/NEG/Invalid for HPIV-4",
            "ROX": "Validity for POS/NEG/Invalid for HPIV-3"            
        }
        
        columns_to_confirm = {
            "FAM": "TRUTHS table for POS/NEG/Invalid for HPIV-1",
            "HEX": "TRUTHS table for POS/NEG/Invalid for HPIV-2",
            "RED647": "TRUTHS table for POS/NEG/Invalid for HPIV-4",
            "ROX": "TRUTHS table for POS/NEG/Invalid for HPIV-3"           
        }
        
        columns_to_pq = {
            "FAM": "PQ Threshold for FAM Rounded RFU Range (HPIV-1)",
            "HEX": "PQ Threshold for HEX Rounded RFU Range (HPIV-2)",
            "RED647": "PQ Threshold for RED647 Rounded RFU Range (HPIV-4)",
            "ROX": "PQ Threshold for ROX Rounded RFU Range (HPIV-3)"
        }

        # Check for invalids
        for chan, validate_column in columns_to_validate.items():
            
            if row[validate_column].lower().find("invalid") > -1:
            
                return "FAIL"
            
        # Check for false positives/negatives
        for chan, confirm_column in columns_to_confirm.items():
            if row[confirm_column].lower().find("false") > -1:
                
                return "FAIL"
            
        # Checks for PQ thresholds
        for chan, pq_value in columns_to_pq.items():
            
            if row[pq_value].lower().find("fail") > -1:
       
                return "FAIL"
            
        return "PASS"

    def stats(self, groupby_settings='Run ID', *args, **kwargs):

        pq_dframe = self.run_pq

        get_stats = {
            "PQ_RESULTS":self.get_pq_results(groupby_settings, pq_dframe),
            "STATS":self.get_pq_results(groupby_settings, pq_dframe)
        }

        return get_stats
        
    def get_pq_results(self, group_stats_by, pq_dframe):
        """ Gets the number of samples that passed PQ and the number that failed """

        all_groups = {}

        pq_categories = {
            'FAM': ['PQ Threshold for FAM Rounded RFU Range (HPIV-1)', 'TRUTHS table for POS/NEG/Invalid for HPIV-1', 'Validity for POS/NEG/Invalid for HPIV-1'],
            'HEX': ['PQ Threshold for HEX Rounded RFU Range (HPIV-2)', 'TRUTHS table for POS/NEG/Invalid for HPIV-2', 'Validity for POS/NEG/Invalid for HPIV-2'],
            'ROX': ['PQ Threshold for ROX Rounded RFU Range (HPIV-3)', 'TRUTHS table for POS/NEG/Invalid for HPIV-3', 'Validity for POS/NEG/Invalid for HPIV-3'],
            'RED647': ['PQ Threshold for RED647 Rounded RFU Range (HPIV-4)', 'TRUTHS table for POS/NEG/Invalid for HPIV-4', 'Validity for POS/NEG/Invalid for HPIV-4']
        }

        for group_cat, series in pq_dframe.groupby([group_stats_by]):

            pq_stats = {
                'POS': {
                    'FAM': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0, 'FALSE_NEG': 0},
                    'HEX': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0, 'FALSE_NEG': 0},
                    'ROX': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0,'FALSE_NEG': 0},
                    'RED647': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0,'FALSE_NEG': 0},
                },
                'NEG':{
                    'FAM': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0, 'FALSE_POS': 0},
                    'HEX': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0, 'FALSE_POS': 0},
                    'ROX': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0, 'FALSE_POS': 0},
                    'RED647': {'COUNT': 0, 'PASS': 0,'FAIL': 0, 'INVALIDS':0, 'FALSE_POS': 0},
                }
            }

            for channel, columns in pq_categories.items():

                data_row_pos = series[series[columns[0]].str.contains('PASS', case=False) & \
                                                              series['SAMPLE Category'].str.contains('POS', case=False) & \
                                                              series['Sample Type'].str.contains('Specimen', case=False) & \
                                                              series[columns[2]].str.contains(r"\bVALID\b")][columns[0]].count()

                data_row_pos_fail = series[series[columns[0]].str.contains('FAIL', case=False) & \
                                                                      series['SAMPLE Category'].str.contains('POS', case=False) & \
                                                                      series['Sample Type'].str.contains('Specimen', case=False) & \
                                                                      series[columns[2]].str.contains(r"\bVALID\b")][columns[0]].count()                                              

                data_row_pos_false = series[series[columns[1]].str.contains('False', case=False) & \
                                                                          series['SAMPLE Category'].str.contains('POS', case=False) & \
                                                                          series['Sample Type'].str.contains('Specimen', case=False) & \
                                                                          series[columns[2]].str.contains(r"\bVALID\b")][columns[1]].count()

                data_row_neg = series[series[columns[0]].str.contains('PASS', case=False) & \
                                                              series['SAMPLE Category'].str.contains('NEG', case=False) & \
                                                              series['Sample Type'].str.contains('Specimen', case=False) & \
                                                              series[columns[2]].str.contains(r"\bVALID\b")][columns[0]].count()

                data_row_neg_fail = series[series[columns[0]].str.contains('FAIL', case=False) & \
                                                                      series['SAMPLE Category'].str.contains('NEG', case=False) & \
                                                                      series['Sample Type'].str.contains('Specimen', case=False) & \
                                                                      series[columns[2]].str.contains(r"\bVALID\b")][columns[0]].count()                                                              

                data_row_neg_false = series[series[columns[1]].str.contains('False', case=False) & \
                                                                          series['SAMPLE Category'].str.contains('NEG', case=False) & \
                                                                          series['Sample Type'].str.contains('Specimen', case=False) & \
                                                                          series[columns[2]].str.contains(r"\bVALID\b")][columns[1]].count()

                data_row_invalid_neg = series[series['SAMPLE Category'].str.contains('NEG', case=False) & \
                                                                             series['Sample Type'].str.contains('Specimen', case=False) & \
                                                                             series[columns[2]].str.contains(r"\bINVALID\b")][columns[0]].count()     

                data_row_invalid_pos = series[series['SAMPLE Category'].str.contains('POS', case=False) & \
                                                                             series['Sample Type'].str.contains('Specimen', case=False) & \
                                                                             series[columns[2]].str.contains(r"\bINVALID\b")][columns[0]].count()

                pq_stats['POS'][channel]['COUNT'] = data_row_pos + data_row_invalid_pos + data_row_pos_fail
                pq_stats['POS'][channel]['PASS'] = data_row_pos
                pq_stats['POS'][channel]['FAIL'] = data_row_pos_fail
                pq_stats['POS'][channel]['INVALIDS'] = data_row_invalid_pos
                pq_stats['POS'][channel]['FALSE_NEG'] = data_row_pos_false

                pq_stats['NEG'][channel]['COUNT'] = data_row_neg + data_row_invalid_neg + data_row_neg_fail
                pq_stats['NEG'][channel]['PASS'] = data_row_neg
                pq_stats['NEG'][channel]['FAIL'] = data_row_neg_fail
                pq_stats['NEG'][channel]['INVALIDS'] = data_row_invalid_neg
                pq_stats['NEG'][channel]['FALSE_POS'] = data_row_neg_false

            all_groups[group_cat] = pq_stats

        return all_groups

    def get_stats_of_valids(self,  group_stats_by, pq_dframe):

        """ Get the stats for 
        (1) VALID (doesn't matter PQ)
        (2) TYPE (pos or neg)

        * ignores positive/negative controls in calculations
        """

        all_group_stats = {}

        categories = {
            'FAM_MEAN':['Validity for POS/NEG/Invalid for HPIV-1', 'FAM-Unrounded RFU Range', 'FAM-Unrounded Ct'],
            'HEX_MEAN':['Validity for POS/NEG/Invalid for HPIV-2', 'HEX-Unrounded RFU Range', 'HEX-Unrounded Ct'],
            'ROX_MEAN':['Validity for POS/NEG/Invalid for HPIV-3', 'ROX-Unrounded RFU Range', 'ROX-Unrounded Ct'],
            'RED647_MEAN':['Validity for POS/NEG/Invalid for HPIV-4', 'RED647-Unrounded RFU Range', 'RED647-Unrounded Ct'],
            'IC_MEAN':['Valid/Invalid for IC', 'IC-Unrounded RFU Range', 'IC-Unrounded Ct']
        }

        for group_cat, series in pq_dframe.groupby([group_stats_by]):

            mean_results = {
                'POS' :{ 
                    'FAM_MEAN': {"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0},
                    'HEX_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0},
                    'ROX_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD": 0},
                    'RED647_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0}
                },
                'NEG' :{
                    'FAM_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0},
                    'HEX_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0},
                    'ROX_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0},
                    'RED647_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0}
                },
                'IC':{
                    'IC_MEAN':{"MEAN_RFU":0, "VALID_COUNT":0, "RFU_STD":0}
                }
            }

            for channel, columns in categories.items():

                if channel == 'IC_MEAN':

                    data_row_ic = series[series[columns[0]].str.contains(r'\bVALID\b', case=False)][columns[1]]

                    ic_mean = data_row_ic.mean() 
                    ic_std =data_row_ic.std() 
                    ic_count = data_row_ic.count()

                    mean_results['IC'][channel]['MEAN_RFU'] = ic_mean
                    mean_results['IC'][channel]['STD_RFU'] = ic_std
                    mean_results['IC'][channel]['VALID_COUNT'] = ic_count

                else:


                    data_row_gen = series[series[columns[0]].str.contains(r'\bVALID\b', case=False) & \
                                                                  series['SAMPLE Category'].str.contains('POS', case=False) & \
                                                                  series['Sample Type'].str.contains('Specimen', case=False)]


                    data_row_rfu_range_pos = data_row_gen[columns[1]].astype(float)
                    data_row_ct_pos = data_row_gen[columns[2]].astype(float).replace("nc", np.NaN)

                    data_row_pos_rfu_mean = data_row_rfu_range_pos.mean()
                    data_row_pos_rfu_count = data_row_rfu_range_pos.count()
                    data_row_pos_rfu_std = data_row_rfu_range_pos.std()

                    data_row_pos_ct_mean = data_row_ct_pos.mean()
                    data_row_pos_ct_count = data_row_ct_pos.count()
                    data_row_pos_ct = data_row_ct_pos.std()

                    data_row_neg_gen = series[series[columns[0]].str.contains(r'\bVALID\b', case=False) &\
                                                                                          series['SAMPLE Category'].str.contains('NEG', case=False) &\
                                                                                          series['Sample Type'].str.contains('Specimen', case=False)]

                    data_row_neg_rfu_range = data_row_neg_gen[[columns[1]]]
                    data_row_neg_ct = data_row_neg_gen[[columns[2]]].replace("nc", np.NaN)

                    data_row_neg_mean_rfu = data_row_neg_rfu_range.mean()
                    data_row_neg_count_rfu = data_row_neg_rfu_range.count()
                    data_row_neg_std_rfu = data_row_neg_rfu_range.std()

                    data_row_neg_mean_ct =  data_row_neg_ct .mean()
                    data_row_neg_std_ct =  data_row_neg_ct .std()

                    mean_results['POS'][channel]['MEAN_RFU'] = data_row_pos_rfu_mean 
                    mean_results['POS'][channel]['VALID_COUNT'] = data_row_pos_rfu_count 
                    mean_results['POS'][channel]['RFU_STD'] = data_row_pos_rfu_std

                    mean_results['POS'][channel]['MEAN_CT'] = data_row_pos_ct_mean 
                    mean_results['POS'][channel]['STD_CT'] =  data_row_pos_ct


                    mean_results['NEG'][channel]['MEAN_RFU'] = data_row_neg_mean_rfu
                    mean_results['NEG'][channel]['VALID_COUNT'] =  data_row_neg_count_rfu
                    mean_results['NEG'][channel]['RFU_STD'] = data_row_neg_std_rfu

                    mean_results['NEG'][channel]['MEAN_CT'] = data_row_neg_mean_ct
                    mean_results['NEG'][channel]['STD_CT'] =  data_row_neg_std_ct

            all_group_stats[group_cat] = mean_results

        return all_group_stats