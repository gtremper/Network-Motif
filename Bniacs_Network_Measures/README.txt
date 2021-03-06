README.txt

BNIACS UPDATED FOLDER:

These are explanations for what each .py file.


___parseNewData.py___ - This parses a ms.pkl, which is multiple sclerosis data.
    
    The ms.pkl's data structure is a dictionary.
    Key: a two-element tuple. 
        First element referrs to the type of matrix that we have as the value.
        The first element is one of the following: 
              ['corr','lcorr,'lacorr','pcorr','pacorr']
        Note: the ms.pkl file we currently have also has
              ['tcorr','tlcorr','tlacorr','tpcorr','tpacorr']
        but does not currently parse for these values.
              
        Second element corresponds to the groups that patients are grouped in.
        In this case, there are four possible groups:
              ["PP","CIS","RR","SP"]
              
    Value: a list of matrices
        for the MS data, the matrices are of size: 70x70
        
        data = load(open("ms.pkl","rb"))
        data[('corr','PP')] returns a list of matrices corresponding to the
                            patient group: 'PP' and corr matrices.
        
        data[('tcorr','PP')] returns a list of matrices that are _thresholded_.
        _thresholded_ means if we do:
            g = nx.DiGraph(data[('tcorr','PP')][0])
            # g will be a directed unweighted graph.
    
    Currently, this parses the current groups we have.
    If we have more types of matrices, add to the list in line 12.
    If we have more group of patients, add to the list in line 13.
    Once you have an updated dictionary, change the FNAME in line 14 to the filename of your dictionary.
    
    main1() computes all the raw measures. (NOTE: Creates a bunch of .pkl files in same directory)
    main2() computes all the t-test and creates T-test charts in the same directory.

___parseHelper.py___ - helper functions for parseNewData.py    
    
___nalz_test.py___ - This makes the t-test charts work. This .py file is used by parseNewData.

    If you want to use this just to make t-test charts, use the: mycompare function.
    
    mycompare - Makes the t-test charts and saves them in same directory.
    
        T-test chart comprises of: Each group compared with one another + each group compared with random-stats.
    
        First arg: is a list of stat measures for each group.
    
            for example, I included in the ms_stats folder a bunch of convertedFormat .pkl files.
            
            d = load(open("ms_stats\ms_convertedFormat_pacorr_undirected.pkl","rb"))
            
            len(d) in this case is 4. (There are 4 patient group types)
            
            Each element in d is of the following data structure:
            [('StronglyConnected','avgoutdeg','stdoutdeg, .... ),([1.0,1.0,1.0,....],[10.0,10.0,...] ...)]
            
            ^ Basically a list of two tuples. 
                First tuple is comprised of all the measure-strings.
                Second tuple is comprised of a list of measure-values correponding to the meausre-strings (index correspondence)
                
            Look in ms_stats\ directory for sample .pkl files.
            
        Second arg: a list of the strings that correspond to the groups in the first argument.
        NOTE: must be in the same order. For example if first arg is 
                [stats for 'PP',stats for 'CIS',stats for 'RR',stats for 'SP']
                Then the second arg must be:
                ['PP','CIS','RR','SP']
                
        Third arg: File for randomized data (randomly generated matrices with stats computed).
            This argument is the same format as an element in d (line 52 of this README).
        
        (look in ms_stats\uprand_D10.pkl and ms_stats\uprand_D10.pkl)
        
        Fourth arg: The name that you want to name the t-test .png file. 
        
        NOTE: this will save in the same directory, so if you have a file with same name, WILL OVERWRITE AND REPLACE.

  GetData(file)
	Takes the name of a pickled dictionary and creats two globally accessible variables: 
	data: the dictionary
	keys: a list of the keys of the dictionary

   RunNetAlgs(dumpData) - computes the graph measures from mynetalgs on the graphs in the appropriate dictionary
	
	To run: First call 'GetData(file)' with the appropriate filename. (Example: aznorbert_corrsd_new.pkl)
	-The file must be a pickled data structure in the form:
	{ "patient ID":   ( {("measure name","corr type"): measure value},  {"matrix type": matrix}, "Patient Type") }

	As in, a dictionary by patient ID, who's values are 3-tuples. First element of which being an originaly empty dictionary,
	second is a dictionary filled with the graphs, both with and with a threshold. And third is the Patient type.
	-"measure names" can be anything as long as they're consistent for each patient
	-"corr type" must be either 'corr', 'lcorr', or 'lacorr'. All 3 for each patient and measure name
	-"matrix type" will inlcude 'corr', 'lcorr', 'lacorr', 'tcorr', 'tlcorr', and 'tlacorr'. The ones starting
	with t mean all edges under a threshold have been removed such that the graph has average
	degree 10 
	-"Patient Type" with be either 'AD', 'MCI', 'NL', or 'CONVERT'

	Argument dumpData specifies whether to dump the described data structure with
	all the computed values. Along with a converted version of the raw data format described
	by the README in Alz_stats. Dumps by defaults

   DegreeNetAlgs(undirected,startdeg,enddeg) 
	Does the same thing as RunNetAlgs but for a range of degree thresholds, starting from
	startdeg (inclusive) to enddeg (inclusive). Defaults to directed graphs of degree 10 only
	Dumps data automatically
	Must run GetData first as with RunNetAlgs. 
   
   RandomDegrees(undirected,startdeg,enddeg)
	Works like DegreeNetAlgs but you don't need to run 'GetData', automatically
	generates 100 random graphs and computes the graph measures on them.
	Dumps data automatically. Defaults to directed graphs of degree 10
        
___mynetalgs.py___ - This file computes MEASURES.
    
    myallmeasures - computes measures of a simple directed (or undirected) graph.
        input: directed or undirected graph.
        
        output: a 2-tuple. 
            First element: ['stronglyConnected','avgoutdeg',...]
            Second element: [.83,1.0,0.0,...]
            
            Both lists are equal in length.
        
___greedy.py___ - Helper functions for mynetalgs' calculation of modularity and linkrank.

___sepGroups.py___ - seperates a PATIENT ID dictionary into the 4 groups NL,AD,MCI, and CONVERTED.

    Structure of Patient ID dictionary:
            Key: PATIENT ID (string)
            
            Value: Three element tuple:
                First element: Dictionary of measures 
                    key: 2 element tuple:
                            first element: String of the measure.
                            second element: String of the type of corr.
                    value: a number.
                
                Second element: Dictionary of matrices
                    key: one of the following - ['corr','lcorr,'lacorr']
                                    or ['tcorr','tlcorr','tlacorr']
                                (t means thresholded)
                    
                    value: a matrix.
                    
                Third element: A string that represents the patient's group.
                    Either "AD" "MCI" "NL" or "CONV".
    
    To seperate this patient ID dictionary into groups: call the sepGroups function with the string
        of your patient ID file name.
        
    Will save as .pkl files in same directory named adCorr,mciCorr,nlCorr,convCorr.
        ^ are basically a list of matrices.
        
___corrData.py___ - fills in the measures for a PATIENT ID dictionary with an empty dictionary of measures.

    main() is the function that does this.
    
    First argument: string of the file name of the PATIENT ID dictionary .pkl file (without the measures)
    
    Second argument: string of the file name that you want the new measure-filled dictionary to me written as.
    
    For example:
    main("aznorbert_corrsd_undirected.pkl","aznorbert_corrsd_new_measures_undirected.pkl")
    
___DataTest.py___ - tests whether two the corrData.py worked. Compares two PATIENT ID dictionaries.
    
    -replace line 10 and 14 with appropiate file names.
    
___graph_helper.py___ - function for graph-swapping. Not used anywhere else.
    

    
    
    
                    
                    
                    
                
        
                    
            



    
    
    
    
    




